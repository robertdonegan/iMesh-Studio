#!/usr/bin/env python3
"""
build_imesh_studio.py
---------------------
Bakes one or more windows of a 1 m LiDAR composite DTM (GeoTIFF, British
National Grid) into a single self-contained HTML file: iMesh Studio.

Usage
-----
    python3 build_imesh_studio.py \
        --tif  path/to/SO83nw_DTM_1m.tif \
        --out  imesh_studio.html \
        --window 382950 235452 2048 "Brook, village and lanes" \
        --window 380500 236452 2048 "Escarpment and field system"

--window takes: EASTING NORTHING SIZE NAME
  EASTING / NORTHING = south-west corner of the window in metres, BNG
  SIZE               = square window side in metres. Use a power of two
                       (1024 / 2048) so the quadtree divides cleanly.

Requires: numpy, tifffile, imagecodecs (for LZW), Pillow.
Elevation is encoded losslessly into the R and G channels of a WebP image
(R = high byte, G = low byte of round((z - zmin) / step)), then embedded as a
data URI. 0.02 m steps keep well inside LiDAR vertical accuracy.
"""

import argparse, base64, datetime, io, json, math, os, re, sys
import numpy as np
import tifffile
from PIL import Image

STEP = 0.02  # metres per encoded unit


# ---------------------------------------------------------------- geodesy
def bng_to_latlon_osgb36(E, N):
    a, b, F0 = 6377563.396, 6356256.909, 0.9996012717
    lat0, lon0, N0, E0 = math.radians(49), math.radians(-2), -100000.0, 400000.0
    e2 = 1 - (b * b) / (a * a)
    n = (a - b) / (a + b)
    lat, M = lat0, 0.0
    for _ in range(40):
        lat = (N - N0 - M) / (a * F0) + lat
        Ma = (1 + n + 1.25 * n * n + 1.25 * n ** 3) * (lat - lat0)
        Mb = (3 * n + 3 * n * n + 2.625 * n ** 3) * math.sin(lat - lat0) * math.cos(lat + lat0)
        Mc = (1.875 * n * n + 1.875 * n ** 3) * math.sin(2 * (lat - lat0)) * math.cos(2 * (lat + lat0))
        Md = (35 / 24) * n ** 3 * math.sin(3 * (lat - lat0)) * math.cos(3 * (lat + lat0))
        M = b * F0 * (Ma - Mb + Mc - Md)
        if abs(N - N0 - M) < 1e-7:
            break
    sl, cl, tl = math.sin(lat), math.cos(lat), math.tan(lat)
    nu = a * F0 / math.sqrt(1 - e2 * sl * sl)
    rho = a * F0 * (1 - e2) / (1 - e2 * sl * sl) ** 1.5
    eta2 = nu / rho - 1
    VII = tl / (2 * rho * nu)
    VIII = tl / (24 * rho * nu ** 3) * (5 + 3 * tl * tl + eta2 - 9 * tl * tl * eta2)
    IX = tl / (720 * rho * nu ** 5) * (61 + 90 * tl * tl + 45 * tl ** 4)
    X = 1 / (cl * nu)
    XI = 1 / (cl * 6 * nu ** 3) * (nu / rho + 2 * tl * tl)
    XII = 1 / (cl * 120 * nu ** 5) * (5 + 28 * tl * tl + 24 * tl ** 4)
    XIIA = 1 / (cl * 5040 * nu ** 7) * (61 + 662 * tl * tl + 1320 * tl ** 4 + 720 * tl ** 6)
    d = E - E0
    return (math.degrees(lat - VII * d ** 2 + VIII * d ** 4 - IX * d ** 6),
            math.degrees(lon0 + X * d - XI * d ** 3 + XII * d ** 5 - XIIA * d ** 7))


def osgb36_to_wgs84(latd, lond):
    a, b = 6377563.396, 6356256.909
    e2 = 1 - (b * b) / (a * a)
    lat, lon = math.radians(latd), math.radians(lond)
    nu = a / math.sqrt(1 - e2 * math.sin(lat) ** 2)
    x = nu * math.cos(lat) * math.cos(lon)
    y = nu * math.cos(lat) * math.sin(lon)
    z = (1 - e2) * nu * math.sin(lat)
    tx, ty, tz, s = 446.448, -125.157, 542.060, 20.4894e-6
    rx, ry, rz = (math.radians(v / 3600) for v in (0.1502, 0.2470, 0.8421))
    X = tx + x * (1 + s) - rz * y + ry * z
    Y = ty + rz * x + y * (1 + s) - rx * z
    Z = tz - ry * x + rx * y + z * (1 + s)
    a2, b2 = 6378137.0, 6356752.3142
    e22 = 1 - (b2 * b2) / (a2 * a2)
    p = math.hypot(X, Y)
    la = math.atan2(Z, p * (1 - e22))
    for _ in range(10):
        nu2 = a2 / math.sqrt(1 - e22 * math.sin(la) ** 2)
        la = math.atan2(Z + e22 * nu2 * math.sin(la), p)
    return math.degrees(la), math.degrees(math.atan2(Y, X))


def bng_to_wgs84(E, N):
    return osgb36_to_wgs84(*bng_to_latlon_osgb36(E, N))


# ---------------------------------------------------------------- raster
def read_geo(tif):
    with tifffile.TiffFile(tif) as t:
        p = t.pages[0]
        tie = p.tags["ModelTiepointTag"].value
        scale = p.tags["ModelPixelScaleTag"].value
    return float(tie[3]), float(tie[4]), float(scale[0]), float(scale[1])


def encode_window(z):
    zmin = float(math.floor(z.min() * 10) / 10)
    v = np.round((z - zmin) / STEP).astype(np.int64)
    if v.max() > 65535:
        raise SystemExit("Elevation range too large for 16-bit encoding — raise STEP.")
    arr = np.dstack([(v >> 8).astype("uint8"), (v & 255).astype("uint8"),
                     np.zeros(v.shape, "uint8")])
    buf = io.BytesIO()
    Image.fromarray(arr).save(buf, format="WEBP", lossless=True, quality=100, method=6)
    raw = buf.getvalue()
    # verify the round trip before we ship it
    back = np.asarray(Image.open(io.BytesIO(raw)).convert("RGB"))
    rt = (back[:, :, 0].astype(np.int32) << 8 | back[:, :, 1].astype(np.int32)) * STEP + zmin
    err = float(np.abs(rt - z).max())
    if err > STEP:
        raise SystemExit(f"WebP round trip lost precision ({err:.4f} m) — switch to PNG.")
    return raw, zmin, err


def slug(s):
    return re.sub(r"[^a-z0-9]+", "-", s.lower()).strip("-")[:40] or "area"


def gridref(E, N):
    """OS 6-figure-ish grid reference for a 100 km square in the SO/SP/ST band."""
    letters = "ABCDEFGHJKLMNOPQRSTUVWXYZ"
    e100, n100 = int(E // 100000), int(N // 100000)
    i500, j500 = (e100 + 10) // 5, (n100 + 5) // 5
    first = letters[(4 - j500) * 5 + i500]
    i, j = (e100 + 10) % 5, (n100 + 5) % 5
    second = letters[(4 - j) * 5 + i]
    return f"{first}{second} {int(E % 100000) // 100:03d}{int(N % 100000) // 100:03d}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tif", required=True)
    ap.add_argument("--template", default=os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                       "imesh_template.html"))
    ap.add_argument("--out", default="imesh_studio.html")
    ap.add_argument("--window", nargs=4, action="append", metavar=("E", "N", "SIZE", "NAME"),
                    required=True)
    ap.add_argument("--note", default="")
    ap.add_argument("--split", metavar="DIR",
                    help="write a Pages-ready folder (index.html + data/*.webp) "
                         "instead of one embedded file")
    args = ap.parse_args()

    print(f"reading {args.tif}")
    dem = tifffile.imread(args.tif).astype(np.float32)
    ox, oy, sx, sy = read_geo(args.tif)
    H, W = dem.shape
    print(f"  {W} x {H} px, origin {ox:.0f} E {oy:.0f} N, {sx:g} m pixels")
    dem[dem < -1e30] = np.nan
    if np.isnan(dem).any():
        print("  filling nodata with nearest finite value")
        med = float(np.nanmedian(dem))
        dem = np.nan_to_num(dem, nan=med)

    regions, total = [], 0
    for E, N, S, name in args.window:
        E, N, S = int(E), int(N), int(S)
        if S & (S - 1):
            print(f"  ! {name}: size {S} is not a power of two; quadtree levels will not divide cleanly")
        c0 = int(round((E - ox) / sx))
        r0 = int(round((oy - (N + S)) / sy))
        if c0 < 0 or r0 < 0 or c0 + S > W or r0 + S > H:
            raise SystemExit(f"window '{name}' ({E},{N},{S}) falls outside the raster "
                             f"({ox:.0f}-{ox+W*sx:.0f} E, {oy-H*sy:.0f}-{oy:.0f} N)")
        z = np.ascontiguousarray(dem[r0:r0 + S, c0:c0 + S])
        raw, zmin, err = encode_window(z)
        total += len(raw)
        lat, lon = bng_to_wgs84(E + S / 2, N + S / 2)
        print(f"  {name}: {S} m square at {gridref(E + S/2, N + S/2)}  "
              f"({lat:.5f}, {lon:.5f})  z {z.min():.2f}–{z.max():.2f} m  "
              f"{len(raw)/1e6:.2f} MB  round-trip error {err*100:.1f} cm")
        regions.append(dict(
            id=slug(name), name=f"{name} — {gridref(E + S/2, N + S/2)}",
            east=E, north=N, w=S, step=STEP, zmin=zmin,
            checkMin=round(float(z.min()), 3), checkMax=round(float(z.max()), 3),
            lat=round(lat, 6), lon=round(lon, 6),
            data="data:image/webp;base64," + base64.b64encode(raw).decode("ascii"),
            _raw=raw,
        ))

    build = dict(
        source=os.path.basename(args.tif),
        date=datetime.date.today().isoformat(),
        note=args.note or ("Environment Agency LiDAR Composite DTM 2022, 1 m, "
                           "Open Government Licence v3."),
    )

    with open(args.template, encoding="utf-8") as f:
        html = f.read()

    if args.split:
        # Pages-ready: small editable index.html, terrain as sibling files
        data_dir = os.path.join(args.split, "data")
        os.makedirs(data_dir, exist_ok=True)
        light = []
        for r in regions:
            fn = f"{r['id']}.webp"
            with open(os.path.join(data_dir, fn), "wb") as f:
                f.write(r.pop("_raw"))
            r.pop("data")
            r["src"] = "data/" + fn
            light.append(r)
        out = os.path.join(args.split, "index.html")
        html = html.replace("/*__REGIONS__*/null", json.dumps(light, indent=1))
        html = html.replace("/*__BUILD__*/null", json.dumps(build, separators=(",", ":")))
        with open(out, "w", encoding="utf-8") as f:
            f.write(html)
        total = sum(os.path.getsize(os.path.join(data_dir, f)) for f in os.listdir(data_dir))
        print(f"\nwrote {out}  ({os.path.getsize(out)/1e3:.0f} kB) "
              f"+ data/ ({total/1e6:.1f} MB, {len(light)} file"
              f"{'s' if len(light) != 1 else ''})")
        print("serve it locally with:  python3 -m http.server 8000")
        return

    for r in regions:
        r.pop("_raw", None)
    html = html.replace("/*__REGIONS__*/null", json.dumps(regions, separators=(",", ":")))
    html = html.replace("/*__BUILD__*/null", json.dumps(build, separators=(",", ":")))
    with open(args.out, "w", encoding="utf-8") as f:
        f.write(html)
    print(f"\nwrote {args.out}  ({os.path.getsize(args.out)/1e6:.1f} MB, "
          f"{len(regions)} area{'s' if len(regions) != 1 else ''})")


if __name__ == "__main__":
    main()
