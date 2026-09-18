#!/usr/bin/env python3
# Figure 3 - Climate of the Eastern Black Sea Region (Türkiye)
# Real data: WorldClim v2.1 monthly climate normals (1970-2000), 5-arc-min.
# Panels: (a) annual precip map, (b) annual temp map, (c) monthly precip,
#         (d) monthly temp, (e) seasonal precip share, + summary.
# Downloads WorldClim (~150 MB) + Natural Earth coast + geoBoundaries provinces.
import requests, zipfile, io, os, numpy as np, warnings
warnings.filterwarnings("ignore")
import rasterio
from rasterio.windows import from_bounds, bounds as wbounds
import geopandas as gpd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MultipleLocator, FuncFormatter

W,S,E,N = 37.2,39.6,42.3,41.75
OUT="climate.png"
BASE="https://geodata.ucdavis.edu/climate/worldclim/2_1/base/"
def fetch_zip(var,res):
    r=requests.get(f"{BASE}wc2.1_{res}_{var}.zip",timeout=280); r.raise_for_status()
    d=f"/tmp/{var}"; os.makedirs(d,exist_ok=True)
    zipfile.ZipFile(io.BytesIO(r.content)).extractall(d); return d
def region_stack(d,var,res):
    arrs=[]; ext=None
    for m in range(1,13):
        with rasterio.open(f"{d}/wc2.1_{res}_{var}_{m:02d}.tif") as ds:
            win=from_bounds(W,S,E,N,ds.transform)
            a=ds.read(1,window=win,masked=True).astype("float32").filled(np.nan)
            if ext is None:
                b=wbounds(win,ds.transform); ext=[b[0],b[2],b[1],b[3]]
            arrs.append(a)
    return np.array(arrs),ext
def sample_pts(d,var,res,pts):
    out={n:[] for n,_,_ in pts}
    for m in range(1,13):
        with rasterio.open(f"{d}/wc2.1_{res}_{var}_{m:02d}.tif") as ds:
            nd=ds.nodata
            for n,lo,la in pts:
                v=float(list(ds.sample([(lo,la)]))[0][0])
                if (nd is not None and v==nd) or (var=="prec" and v<0) or v<-1e29: v=np.nan
                out[n].append(v)
    return out

pts=[("Rize",40.52,41.02),("Trabzon",39.72,41.00),("Giresun",38.39,40.91),
     ("Artvin",41.82,41.18),("Gumushane",39.48,40.46)]
dp=fetch_zip("prec","5m"); dt=fetch_zip("tavg","5m")
P,extP=region_stack(dp,"prec","5m"); T,extT=region_stack(dt,"tavg","5m")
ann_p=np.nansum(P,axis=0); ann_p[np.all(np.isnan(P),axis=0)]=np.nan
ann_t=np.nanmean(T,axis=0)
pst=sample_pts(dp,"prec","5m",pts); tst=sample_pts(dt,"tavg","5m",pts)
reg_month=np.nanmean(P.reshape(12,-1),axis=1)
seas={"Winter (DJF)":reg_month[[11,0,1]].sum(),"Spring (MAM)":reg_month[2:5].sum(),
      "Summer (JJA)":reg_month[5:8].sum(),"Autumn (SON)":reg_month[8:11].sum()}
tot=sum(seas.values()); seas_pct={k:float(v/tot*100) for k,v in seas.items()}
mon=["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"]

def gj(url):
    r=requests.get(url,timeout=120); open("/tmp/x.geojson","wb").write(r.content)
    return gpd.read_file("/tmp/x.geojson",bbox=(W,S,E,N))
coast=gj("https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/geojson/ne_10m_coastline.geojson")
prov=gj("https://media.githubusercontent.com/media/wmgeolab/geoBoundaries/main/releaseData/gbOpen/TUR/ADM1/geoBoundaries-TUR-ADM1.geojson")

plt.rcParams.update({"font.family":"DejaVu Sans","font.size":10})
asp=1/np.cos(np.radians(0.5*(S+N)))
fig=plt.figure(figsize=(16,9.6))
gs=GridSpec(2,3,height_ratios=[1.15,1],hspace=0.30,wspace=0.24)
def basemap(ax):
    coast.plot(ax=ax,color="black",linewidth=0.7,zorder=5)
    prov.boundary.plot(ax=ax,color="0.35",linewidth=0.4,linestyle="--",alpha=0.6,zorder=4)
    for n,lo,la in pts:
        ax.plot(lo,la,"o",ms=4,mfc="black",mec="white",mew=0.6,zorder=6)
        ax.annotate(n,(lo,la),xytext=(0,4),textcoords="offset points",ha="center",fontsize=7.5,zorder=7)
    ax.set_xlim(W,E); ax.set_ylim(S,N); ax.set_aspect(asp)
    ax.xaxis.set_major_locator(MultipleLocator(1)); ax.yaxis.set_major_locator(MultipleLocator(1))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v,_:f"{v:g}°E"))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v,_:f"{v:g}°N")); ax.tick_params(labelsize=8)
axa=fig.add_subplot(gs[0,0]); axa.set_facecolor("#DCEAF4")
plev=[400,600,800,1000,1200,1400,1600,2000,2500,3000]; cmapP=plt.get_cmap("YlGnBu"); normP=BoundaryNorm(plev,cmapP.N,extend="both")
imP=axa.imshow(ann_p,extent=extP,origin="upper",cmap=cmapP,norm=normP,zorder=1); basemap(axa)
axa.set_title("(a) Mean annual precipitation",fontsize=11,loc="left")
cbP=fig.colorbar(imP,ax=axa,fraction=0.046,pad=0.02,extend="both"); cbP.set_label("mm yr$^{-1}$",fontsize=8); cbP.ax.tick_params(labelsize=7)
axb=fig.add_subplot(gs[0,1]); axb.set_facecolor("#DCEAF4")
tlev=np.arange(-2,20.1,2); cmapT=plt.get_cmap("RdYlBu_r"); normT=BoundaryNorm(tlev,cmapT.N,extend="both")
imT=axb.imshow(ann_t,extent=extT,origin="upper",cmap=cmapT,norm=normT,zorder=1); basemap(axb)
axb.set_title("(b) Mean annual temperature",fontsize=11,loc="left")
cbT=fig.colorbar(imT,ax=axb,fraction=0.046,pad=0.02,extend="both"); cbT.set_label("°C",fontsize=8); cbT.ax.tick_params(labelsize=7)
axs=fig.add_subplot(gs[0,2]); axs.axis("off")
txt=("Climate summary (WorldClim v2.1)\n\n"
     f"• Annual precipitation: {np.nanmin(ann_p):.0f}-{np.nanmax(ann_p):.0f} mm yr$^{{-1}}$\n"
     "• Wettest on the coastal ranges, driest inland\n"
     f"• Mean annual temperature: {np.nanmin(ann_t):.0f}-{np.nanmax(ann_t):.0f} °C\n"
     f"• Coastal regime: autumn-early-winter maximum\n  (Rize wettest in {mon[int(np.nanargmax(pst['Rize']))]})\n\n"
     "Normals 1970-2000, ~9 km grid.\nStations sampled at provincial centres.\nProjection: geographic, WGS84.")
axs.text(0.02,0.98,txt,va="top",ha="left",fontsize=9,linespacing=1.5,bbox=dict(boxstyle="round,pad=0.6",fc="#F5F5F0",ec="0.6"))
mcol={"Rize":"#1f77b4","Trabzon":"#2ca02c","Giresun":"#ff7f0e","Artvin":"#9467bd","Gumushane":"#8c564b"}
months=range(1,13); ml=[m[:3] for m in mon]
axc=fig.add_subplot(gs[1,0])
for n,_,_ in pts: axc.plot(months,pst[n],"-o",ms=3,lw=1.3,color=mcol[n],label=n)
axc.set_title("(c) Monthly precipitation",fontsize=11,loc="left"); axc.set_ylabel("mm")
axc.set_xticks(list(months)); axc.set_xticklabels(ml,fontsize=7,rotation=45); axc.grid(lw=0.4,color="0.85")
axc.legend(fontsize=7,ncol=2,frameon=False); axc.yaxis.set_minor_locator(MultipleLocator(25))
axd=fig.add_subplot(gs[1,1])
for n,_,_ in pts: axd.plot(months,tst[n],"-o",ms=3,lw=1.3,color=mcol[n],label=n)
axd.set_title("(d) Monthly mean temperature",fontsize=11,loc="left"); axd.set_ylabel("°C")
axd.set_xticks(list(months)); axd.set_xticklabels(ml,fontsize=7,rotation=45); axd.grid(lw=0.4,color="0.85"); axd.axhline(0,color="0.5",lw=0.6)
axe=fig.add_subplot(gs[1,2])
labels=list(seas_pct.keys()); vals=[seas_pct[k] for k in labels]; pcol=["#4575b4","#91cf60","#fc8d59","#8073ac"]
axe.pie(vals,labels=[f"{k}\n{v:.0f}%" for k,v in zip(labels,vals)],colors=pcol,startangle=90,counterclock=False,
        textprops={"fontsize":8},wedgeprops={"ec":"white","lw":1}); axe.set_title("(e) Seasonal precipitation share",fontsize=11,loc="left")
fig.suptitle("Climate of the Eastern Black Sea Region (Türkiye)",fontsize=14,fontweight="bold",y=0.98)
fig.text(0.5,0.012,"Data: WorldClim v2.1 climate normals (1970-2000). Coast & provinces: Natural Earth / geoBoundaries. WGS84.",ha="center",fontsize=8,style="italic",color="0.4")
fig.savefig(OUT,dpi=300,bbox_inches="tight",facecolor="white")
print("saved",OUT)
