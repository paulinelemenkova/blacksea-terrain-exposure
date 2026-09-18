# Input data

All datasets used by this project are **open-access**. Most are retrieved
automatically at run time; two vector inputs and the geological map must be
placed here manually because they are redistributed under their own licences.

## Retrieved automatically by the scripts
- **Relief / bathymetry** — SRTM/GEBCO composite, via the GMT data server
  (`@earth_relief_03s`).
- **Climate normals** — WorldClim v2.1 (downloaded by `climate.py`).
- **Forest change** — Hansen Global Forest Change v1.11 tiles
  (downloaded by `forest_dynamics.sh`).
- **Population** — WorldPop 2020, ~1 km (downloaded by `exposure.sh`).
- **Coast / rivers** — GSHHG, via GMT.

## Place here manually
- `data/geology/` — the 1:500,000 Geological Map of Turkiye
  (MTA, 1961; digital vector version digitised by Orr, 2002), as MapInfo TAB:
  `Turkey 500k Lithology V1.TAB`, `Turkey 500k Faults V1.TAB`,
  `Turkey 500k MinOcc V1.TAB`. Used by `geology_maps.py`.
- `data/ne_coast.geojson`, `data/ne_adm1.geojson` — Natural Earth 1:10m
  coastline and admin-1 boundaries (used by `geology_maps.py`).

Nothing in this folder is tracked by git (see `.gitignore`).
