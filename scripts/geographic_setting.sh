#!/usr/bin/env bash
# =============================================================================
# Figure 1 - Geographic setting, Eastern Black Sea Region (Turkiye)  [rev.2]
# Fixes: 3" relief (sharper) · no white strip over inset · repositioned
#        aliceblue non-bold sea label · ASCII-safe names · white graticule · WEsN
# Requires GMT 6.x + GDAL. Run:  bash geographic_setting.sh
# =============================================================================
set -e
REG=-R37.2/42.3/39.6/41.75
PROJ=-JM17c

# --- relief: 3 arc-sec (~90 m) for sharpness. For your local global GEBCO use:
# gmt grdcut /Volumes/TOSHIBA/DATA/GEBCO_Marmara/GEBCO_2026.nc $REG -Grelief.nc
gmt grdcut @earth_relief_03s $REG -Grelief.nc

# --- province boundaries: local file if present, else geoBoundaries TUR ADM1
if [ -f geoBoundaries-TUR-ADM1.gmt ]; then cp geoBoundaries-TUR-ADM1.gmt tur.gmt
else
  curl -sL --max-time 90 -o tur.geojson \
    https://media.githubusercontent.com/media/wmgeolab/geoBoundaries/main/releaseData/gbOpen/TUR/ADM1/geoBoundaries-TUR-ADM1.geojson
  ogr2ogr -f OGR_GMT tur.gmt tur.geojson
fi

gmt makecpt -Cgeo -T-2500/3700 > topo.cpt

gmt begin geographic_setting png E500
  gmt set PS_CHAR_ENCODING ISOLatin1+
  gmt set FONT_ANNOT_PRIMARY 10p,Helvetica,black FONT_LABEL 11p,Helvetica,black
  gmt set MAP_FRAME_TYPE fancy MAP_FRAME_WIDTH 3.5p MAP_FRAME_PEN 1p
  gmt set FORMAT_GEO_MAP dddF MAP_GRID_PEN_PRIMARY 0.4p,white   # thin white graticule

  # relief + white graticule (g0.5); annotate W E N, tick s  (WEsN)
  gmt grdimage relief.nc -Ctopo.cpt -I+d $REG $PROJ -Bxa1f0.5g0.5 -Bya1f0.5g0.5 -BWEsN

  gmt coast $REG $PROJ -Df -W0.5p,black -N1/0.9p,gray20 -I1/0.8p,steelblue -I2/0.4p,steelblue
  gmt plot tur.gmt -W0.8p,red1

  gmt plot -Sc0.14c -Gred3 -W0.5p,black <<'EOF'
37.88 40.98
38.39 40.91
39.72 41.00
40.52 41.02
41.82 41.18
39.48 40.46
40.23 40.26
EOF
  # ASCII-safe names (Turkish s-cedilla / dotted-I are absent from the PS encoding)
  gmt text -F+f9.5p,Helvetica-Bold,black+jTC -Gwhite@30 -D0/-0.16c <<'EOF'
37.88 40.98 Ordu
38.39 40.91 Giresun
39.72 41.00 Trabzon
40.52 41.02 Rize
41.82 41.18 Artvin
39.48 40.46 Gumushane
40.23 40.26 Bayburt
EOF
  gmt text -F+f8.5p,Helvetica-BoldOblique,khaki+jCM <<'EOF'
38.441 40.566 GIRESUN
39.781 40.809 TRABZON
40.816 40.920 RIZE
41.848 41.039 ARTVIN
EOF

  # sea label: right of inset, smaller, non-bold, aliceblue
  echo "40.2 41.58 BLACK SEA" | gmt text -F+f14p,Helvetica-Oblique,aliceblue+jCM

  gmt basemap -LjBL+w100k+o0.6c/0.6c+f+u --FONT_ANNOT_PRIMARY=9p
  gmt basemap -TdjTR+w1.1c+o0.5c/0.5c+l,,,N
  gmt colorbar -Ctopo.cpt -DJBC+o0/1.1c+w12c/0.35c+h -Bxa1000f500+l"Elevation / bathymetry (m)"

  # locator inset: sea-filled panel removes the white strip
  gmt inset begin -DjTL+w5c/2.4c+o0.3c/0.3c -F+glightblue+p0.6p,black
    gmt coast -R25.5/45/35.5/42.6 -JM5c -Ggray85 -Slightblue -N1/0.5p,gray30 -Wfaint -A2000
    echo "37.2 39.6 42.3 41.75" | gmt plot -Sr+s -W1.2p,red3
    echo "29.5 39.2 TÜRKIYE" | gmt text -F+f8p,Helvetica-Bold,gray20+jLM -N
  gmt inset end

  gmt text -R0/1/0/1 -JX17c/11c -F+f7.5p,Helvetica-Oblique,gray30+jBR -N -Gwhite@40 <<'EOF'
1.0 -0.075 Relief: SRTM15+/GEBCO 3" (GMT). Coast & rivers: GSHHG. Provinces: geoBoundaries TUR ADM1. Proj: Mercator, WGS84.
EOF
gmt end show
