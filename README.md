# iMesh Studio

Quadtree mesh-refinement plates from real LiDAR terrain, for the Flood Modeller
v8 iMesh keynote. Dial a refinement slider, export frame sequences for
After Effects.

## Run it

```
python3 -m http.server 8000
# then open http://localhost:8000
```

A server is required — the tool reads elevation out of a canvas with
`getImageData`, and Chrome taints the canvas for images loaded over `file://`.
Opening `index.html` directly will show a "terrain failed to load" message.

## Layout

```
index.html                  the whole tool — this is the file you edit
data/*.webp                 terrain, one file per area (elevation in R+G)
tools/build_imesh_studio.py cuts new areas from a LiDAR GeoTIFF
tools/imesh_template.html   index.html with no terrain baked in
```

`index.html` and `tools/imesh_template.html` are the same file apart from one
JSON blob near the top of the script. **Edit `index.html`.** When you want to
change which areas ship, copy your edits back into the template and re-run the
build — or just hand-edit the `REGIONS` array in `index.html`, which is the
first thing in the `<script>`.

## Adding an area

Download a 1 m composite DTM tile from the DEFRA Survey Data Download service,
then:

```
pip install numpy tifffile imagecodecs pillow

python3 tools/build_imesh_studio.py \
  --tif SO84se_DTM_1m.tif \
  --template tools/imesh_template.html \
  --split . \
  --window 385000 240000 2048 "Upton upon Severn"
```

`--window` takes `EASTING NORTHING SIZE NAME` — south-west corner in metres
BNG, square side in metres. Use a power of two. Repeat `--window` for more
areas. Each 2048 m square adds about 1.3 MB to `data/`.

Drop `--split .` and use `--out imesh_studio.html` instead to get a single
self-contained file that runs from `file://` with no server. Useful for sending
to someone who will not run a terminal.

## How refinement works

A cell subdivides when the elevation range inside it, multiplied by a class
weight, exceeds a tolerance. Cells nest 16 → 8 → 4 → 2 → 1 m.

Weights come from the DTM itself: slope and elevation change, built-up areas
(detected from the density of kerb- and wall-scale structure), and field
texture. River corridors get a **hard cap** rather than a weight — channels
cannot go finer than the chosen size wherever the slider sits. Weighting alone
was not enough, because channel banks have real relief and kept refining.

The slider is calibrated in mesh density, not tolerance: equal slider steps give
equal steps in log cell count.

## Exporting for After Effects

Set **Output → Background: Mesh only (alpha)** to composite over your own
basemap. From 0 → To 0.85 over 90 frames with easing is a good default — above
about 0.92 the last uncapped cells all flip to 1 m at once and the end of the
sweep looks abrupt.

The exported zip contains `frames.txt` mapping each frame to its refinement
value, so a faked slider in AE can key against the real numbers. **Cell counts
.csv** gives cells and % reduction across 101 steps for captions.

## Basemaps

Hillshade and Ground elevation are generated from the DTM in-browser and always
work. Satellite (Esri), OSM and Carto pull web tiles; if a tile server refuses
CORS the canvas is tainted and PNG export is blocked — the tool says so. Esri
World Imagery is the most reliable of the three.

## Attribution

Terrain: Environment Agency LiDAR Composite DTM 2022, 1 m. Contains public
sector information licensed under the Open Government Licence v3.
