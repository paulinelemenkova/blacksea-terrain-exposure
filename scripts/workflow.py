#!/usr/bin/env python3
# Figure 7 (workflow) - Reproducible open-data conditioning-factor workflow (Path B).
# No landslide inventory, no ML model, no SHAP/hazard/hotspot/scenario, no validation.
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle, FancyArrowPatch
from matplotlib.transforms import Bbox
plt.rcParams.update({"font.family":"sans-serif",
  "font.sans-serif":["Nimbus Sans","Arial","Liberation Sans","Arimo","Helvetica","DejaVu Sans"]})

fig=plt.figure(figsize=(15,8.6)); ax=fig.add_axes([0,0,1,1]); ax.axis("off")
ax.set_xlim(0,1); ax.set_ylim(0,1)
ax.add_patch(Rectangle((0,0.38),1,0.585,fc="#f7f8fa",ec="none",zorder=0))
BS="round,pad=0.004,rounding_size=0.012"

def card(x,y,w,h,color,title,items,ncol=1,tsize=10.5,isize=7.6,istep=0.030):
    ax.add_patch(FancyBboxPatch((x+0.003,y-0.005),w,h,boxstyle=BS,fc="#0000000e",ec="none",zorder=1))
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle=BS,fc="white",ec="#d3d7dd",lw=1.0,zorder=2))
    ax.add_patch(Rectangle((x+0.004,y+0.010),0.0055,h-0.020,fc=color,ec="none",zorder=3))
    ax.text(x+w/2+0.006,y+h-0.020,title,fontsize=tsize,fontweight="bold",color=color,va="top",ha="center",zorder=4)
    ax.plot([x+0.022,x+w-0.012],[y+h-0.050,y+h-0.050],color=color,lw=1.0,alpha=0.45,zorder=4)
    per=int(np.ceil(len(items)/ncol)); colw=(w-0.028)/ncol
    for ci in range(ncol):
        sub=items[ci*per:(ci+1)*per]; ty=y+h-0.066; tx=x+0.018+ci*colw
        for it in sub:
            ax.text(tx,ty,"•  "+it,fontsize=isize,color="#2b2f36",va="top",ha="left",zorder=4); ty-=istep

def arrow(x0,y0,x1,y1,color="#8a9099",lw=2.4):
    ax.add_patch(FancyArrowPatch((x0,y0),(x1,y1),arrowstyle="-|>",mutation_scale=16,lw=lw,color=color,zorder=1,shrinkA=0,shrinkB=0))

COL=["#2166ac","#4393c3","#1b9e77","#5aa02c","#7570b3","#d95f02"]
xs=[0.018,0.183,0.348,0.513,0.678,0.845]; w=0.150; ytop=0.60; htop=0.26
titles=["1 · Data sources","2 · Pre-processing","3 · Geomorphometry","4 · Terrain analysis","5 · Conditioning factors","6 · Outputs & release"]
items=[
 ["DEM (SRTM / GEBCO)","Climate (WorldClim v2.1)","Geology, faults (MTA 1:500k)","Forest (Hansen GFC v1.11)","Coast, boundaries (Nat. Earth)","(open access; no inventory)"],
 ["Reproject to WGS84 grid","Clip to study area","Resample to common grid","Depression filling (routing)","Layer stacking"],
 ["Elevation, slope, aspect","Plan & profile curvature","Terrain Ruggedness Index","Topographic Wetness Index","Hypsometry & relief stats"],
 ["Stream Power Index (SPI)","Sediment Transport Index","Topographic Position Index","Weiss landform classes","Drainage network, profile"],
 ["Lithological–structural frame","Climatic regime (rain, temp.)","Forest dynamics 2000–2020","Distance to faults / rivers","Synthesis of predisposing controls"],
 ["Characterisation maps & figures","Open geospatial data baseline","Documented GMT/Python scripts","Inventory-independent result","Basis for future susceptibility"],
]
for i in range(6):
    card(xs[i],ytop,w,htop,COL[i],titles[i],items[i])
    if i<5: arrow(xs[i]+w+0.002,ytop+htop/2,xs[i+1]-0.002,ytop+htop/2)

# foundation band: reproducibility & tools (underpins all stages; not a pipeline step)
by,bh=0.40,0.15; bx,bw=0.018,0.977
card(bx,by,bw,bh,"#455a64","Reproducibility & tools  (underpin all stages)",
     ["Open-access inputs only — no proprietary data, no landslide inventory",
      "GMT 6 for relief retrieval and mapping",
      "Python: numpy, scipy, rasterio, geopandas, matplotlib",
      "GDAL for raster / vector I/O and reprojection",
      "Deterministic, end-to-end reproducible from public sources",
      "All scripts released with the paper"],
     ncol=3, tsize=11, isize=8.4, istep=0.052)
# light connectors from the pipeline down to the foundation
for i in (0,5): arrow(xs[i]+w/2, ytop-0.004, xs[i]+w/2, by+bh+0.006, color="#90a4ae", lw=1.6)

fig.text(0.5,0.94,"Reproducible open-data workflow — terrain and conditioning-factor characterisation, Eastern Black Sea Region",
         ha="center",fontsize=13.5,fontweight="bold",color="#1a1d22")
fig.text(0.5,0.912,"From open-access data through geomorphometric and terrain analysis to a released, inventory-independent conditioning-factor baseline",
         ha="center",fontsize=9.5,color="#5b6069",style="italic")
fig.savefig("figures/workflow.png",dpi=200,bbox_inches=Bbox([[0,3.3],[15,8.35]]),facecolor="#f7f8fa")
print("saved")
