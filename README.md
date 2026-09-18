# Eastern Black Sea — terrain conditioning factors and population exposure

An open, fully reproducible **Python + GMT** workflow that characterises the
geomorphometric and environmental conditioning factors of slope instability in
the Eastern Black Sea Region of Türkiye, and links them to **population
exposure**. Every layer and figure is generated from open-access data by the
scripts in this repository — no proprietary data and no landslide inventory are
used. The heuristic predisposition index produced here is an inventory-independent
screening product, **not** a validated susceptibility or hazard map.

## What it does

From open relief, climate, geology, forest and population data the workflow derives:

- **Geomorphometry** — elevation, slope, nine-class aspect, plan/profile curvature,
  Terrain Ruggedness Index, Topographic Wetness Index, hypsometry.
- **Hydro-geomorphometry** — Stream Power and Sediment Transport indices,
  Topographic Position Index, six-class Weiss landforms, drainage network, swath profile.
- **Lithological–structural framework** — rock-type associations, faults, mineral occurrences.
- **Climatic regime** — precipitation/temperature normals and seasonality.
- **Forest dynamics** — two decades of Hansen Global Forest Change.
- **Predisposition synthesis + exposure** — a weighted-overlay predisposition index
  intersected with gridded population to quantify who lives on the most predisposed terrain.

## Repository layout

```
blacksea-terrain-exposure/
├── README.md
├── LICENSE                 # MIT
├── CITATION.cff
├── requirements.txt        # Python dependencies
├── environment.yml         # conda env (includes GMT)
├── .gitignore
├── data/                   # input data (not tracked) — see data/README.md
└── scripts/                # analysis + figure scripts
```

## Scripts

Each script is self-contained and writes its output to `figures/`.
`.sh` scripts drive the Generic Mapping Tools (GMT) directly and embed Python for
raster analysis; `.py` scripts use the Python geospatial stack. Every script's
header documents its specific inputs.

| Script | Output | Main data |
|---|---|---|
| `geographic_setting.sh` | `geographic_setting.png` | SRTM/GEBCO, GSHHG, provinces |
| `geology_maps.py`       | `geological_setting.png`, `lithological_setting.png` | MTA 1:500,000 (Orr 2002) |
| `climate.py`            | `climate.png` | WorldClim v2.1 |
| `geomorphology.sh`      | `geomorphology.png` | SRTM/GEBCO |
| `terrain_analysis.sh`   | `terrain_analysis.png` | SRTM/GEBCO |
| `forest_dynamics.sh`    | `forest_dynamics.png` | Hansen GFC v1.11 |
| `synthesis.py`          | `synthesis.png` | derived conditioning factors |
| `exposure.sh`           | `exposure.png` | derived index + WorldPop 2020 |
| `workflow.py`           | `workflow.png` | schematic (no external data) |

## Installation

**conda (recommended, installs GMT too):**

```bash
conda env create -f environment.yml
conda activate blacksea-terrain
```

**pip (Python stack only; install GMT >= 6.4 separately):**

```bash
pip install -r requirements.txt
```

Core libraries: `numpy`, `scipy`, `rasterio`, `geopandas`, `shapely`, `pyproj`,
`pyogrio`, `matplotlib`, `requests`, `pillow`; plus **GMT >= 6.4** for the `.sh` scripts.

## Usage

Run any script from the repository root; figures are written to `figures/`.

```bash
mkdir -p figures
python scripts/climate.py
bash   scripts/exposure.sh
```

Most inputs are fetched automatically at run time (SRTM/GEBCO via GMT, WorldClim,
Hansen GFC, WorldPop). The geological map and the Natural Earth vectors must be
placed under `data/` first — see [`data/README.md`](data/README.md).

## Data availability

All inputs are open-access and are documented in `data/README.md`. This repository
is under active development; a tagged release will be archived on Zenodo with a
permanent DOI.

## Reproducibility notes

- All layers share a common WGS84 grid (~275 m) clipped to the study window
  (37.2–42.3° E, 39.6–41.75° N); land and sea are separated by the sign of elevation.
- Flow routing uses a heap-based priority-flood fill followed by D8 accumulation,
  so drainage is continuous to the coast.
- The predisposition index is a transparent, expert-weighted overlay
  (slope 0.35, wetness 0.25, plan-curvature convergence 0.20, stream power 0.20),
  classified into five quantiles; it is not calibrated against observed failures.

## License

Released under the MIT License (see `LICENSE`).
