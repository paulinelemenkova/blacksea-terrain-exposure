#!/usr/bin/env bash
# Figure 8 - Two decades of forest dynamics, Eastern Black Sea Region (Türkiye), 2000-2020
# REAL data: Hansen/UMD Global Forest Change v1.11 (30 m). Right column: compact 2x3 legend
# (no title) + (G) forest-cover change by elevation band (loss/gain km2 + net line, REAL
# Hansen numbers); stats & findings moved to bottom annotation.
set -e
REG=-R37.2/42.3/39.6/41.75
gmt grdcut @earth_relief_03s $REG -Gdem03.nc
gmt grdsample dem03.nc -I0.0025 -Gdem.nc
gmt grdconvert dem.nc dem.tif=gd:GTiff
curl -sL -o coast.geojson https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/geojson/ne_10m_coastline.geojson
python3 << 'PY'
import os,warnings; warnings.filterwarnings("ignore")
os.environ["GDAL_DISABLE_READDIR_ON_OPEN"]="EMPTY_DIR"; os.environ["CPL_VSIL_CURL_USE_HEAD"]="NO"
import numpy as np, rasterio, geopandas as gpd
from rasterio.warp import reproject, Resampling
from scipy.ndimage import binary_dilation
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator, FuncFormatter, AutoMinorLocator
plt.rcParams.update({"font.family":"sans-serif",
  "font.sans-serif":["Nimbus Sans","Arial","Liberation Sans","Arimo","Helvetica","DejaVu Sans"],"font.size":9})
W,E,S,N=37.2,42.3,39.6,41.75
dd=rasterio.open("dem.tif"); dem=dd.read(1).astype("float32"); dem[dem<-1e4]=np.nan
ny,nx=dem.shape; tr=dd.transform; crs=dd.crs; ext=[W,E,S,N]
latc=0.5*(S+N); cell=(dd.res[0]*111320*np.cos(np.radians(latc)))*(dd.res[1]*110540)/1e6; land=dem>=0
base="https://storage.googleapis.com/earthenginepartners-hansen/GFC-2023-v1.11/"
tb={"50N_030E":(30,40,40,50),"50N_040E":(40,40,50,50),"40N_030E":(30,30,40,40),"40N_040E":(40,30,50,40)}
TH=30
cf=np.zeros((ny,nx),'f4'); cl20=np.zeros((ny,nx),'f4'); cl10=np.zeros((ny,nx),'f4'); cg=np.zeros((ny,nx),'f4')
f2000=lo20=lo10=gain=0.0
for t,(l,b,r,tp) in tb.items():
    bb=(max(W,l),max(S,b),min(E,r),min(N,tp))
    if bb[0]>=bb[2] or bb[1]>=bb[3]: continue
    latm=0.5*(bb[1]+bb[3]); nat=(0.00025*111320*np.cos(np.radians(latm)))*(0.00025*110540)/1e6
    dtc=rasterio.open("/vsicurl/"+base+f"Hansen_GFC-2023-v1.11_treecover2000_{t}.tif")
    dly=rasterio.open("/vsicurl/"+base+f"Hansen_GFC-2023-v1.11_lossyear_{t}.tif")
    dgn=rasterio.open("/vsicurl/"+base+f"Hansen_GFC-2023-v1.11_gain_{t}.tif")
    w=dtc.window(*bb); wt=dtc.window_transform(w)
    tc=dtc.read(1,window=w); ly=dly.read(1,window=dly.window(*bb)); gn=dgn.read(1,window=dgn.window(*bb))
    n=min(tc.shape[0],ly.shape[0],gn.shape[0]); m=min(tc.shape[1],ly.shape[1],gn.shape[1]); tc,ly,gn=tc[:n,:m],ly[:n,:m],gn[:n,:m]
    f=(tc>=TH); l20=f&(ly>=1)&(ly<=20); l10=f&(ly>=1)&(ly<=10); g=(gn==1)&~f
    f2000+=f.sum()*nat; lo20+=l20.sum()*nat; lo10+=l10.sum()*nat; gain+=g.sum()*nat
    for src,acc in [(f,cf),(l20,cl20),(l10,cl10),(g,cg)]:
        tmp=np.zeros((ny,nx),'f4'); reproject(src.astype('f4'),tmp,src_transform=wt,src_crs=dtc.crs,dst_transform=tr,dst_crs=crs,resampling=Resampling.average,dst_nodata=0.0); acc+=tmp
    for d in (dtc,dly,dgn): d.close()
cf=np.clip(cf,0,1); cl20=np.clip(cl20,0,1); cl10=np.clip(cl10,0,1); cg=np.clip(cg,0,1)
a00=f2000; a10=f2000-lo10; a20=f2000-lo20+gain; net=gain-lo20; pnet=100*net/a00
bands=[(0,500),(500,1000),(1000,1500),(1500,2000),(2000,2500),(2500,9000)]; blab=["0–\n500","500–\n1000","1000–\n1500","1500–\n2000","2000–\n2500",">2500"]
decl=[]; loss_b=[]; gain_b=[]; net_b=[]
for lo,hi in bands:
    bm=land&(dem>=lo)&(dem<hi); f0=float((cf*bm).sum())*cell
    lb=float((cl20*bm).sum())*cell; gb=float((cg*bm).sum())*cell
    decl.append(100*lb/f0 if f0>0 else 0.0); loss_b.append(lb); gain_b.append(gb); net_b.append(gb-lb)
coast=gpd.read_file("coast.geojson",bbox=(W,S,E,N)); AR=1/np.cos(np.radians(latc))
CIT=[("Giresun",38.39,40.91),("Trabzon",39.72,41.00),("Rize",40.52,41.02),("Artvin",41.82,41.18)]
FG="#2e7d32"; NF="#efe8c8"; STB="#bfe0a8"; LS="#e31a1c"; GN="#1f78b4"
fig=plt.figure(figsize=(16.5,8.2))
gs=GridSpec(2,4,figure=fig,width_ratios=[1,1,1,0.80],height_ratios=[1,1],left=0.035,right=0.985,top=0.935,bottom=0.20,hspace=0.17,wspace=0.22)
def basemap(ax):
    coast.plot(ax=ax,color="black",lw=0.5,zorder=5); ax.set_xlim(W,E); ax.set_ylim(S,N); ax.set_aspect(AR)
    for xln in range(38,43): ax.axvline(xln,color="white",lw=0.4,alpha=0.75,zorder=4)
    for yln in range(40,42): ax.axhline(yln,color="white",lw=0.4,alpha=0.75,zorder=4)
    for n,lo,la in CIT: ax.plot(lo,la,"o",ms=3.2,mfc="black",mec="white",mew=0.5,zorder=6)
    ax.xaxis.set_major_locator(MultipleLocator(1)); ax.yaxis.set_major_locator(MultipleLocator(1))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x,_:f"{x:g}°E")); ax.yaxis.set_major_formatter(FuncFormatter(lambda x,_:f"{x:g}°N")); ax.tick_params(labelsize=8)
def fmap(ax,cover,title):
    ax.set_facecolor("#DCEAF4"); ax.imshow(np.where(~land,np.nan,(cover>=0.5).astype(float)),extent=ext,origin="upper",cmap=ListedColormap([NF,FG]),vmin=0,vmax=1); basemap(ax)
    ax.set_title(title,fontsize=11,loc="left",fontweight="bold")
axA=fig.add_subplot(gs[0,0]); fmap(axA,cf,"(A) Forest extent – 2000")
axB=fig.add_subplot(gs[0,1]); fmap(axB,cf-cl10,"(B) Forest extent – 2010")
axC=fig.add_subplot(gs[0,2]); fmap(axC,cf-cl20,"(C) Forest extent – 2020")
axD=fig.add_subplot(gs[1,0]); axD.set_facecolor("#DCEAF4")
lossm=binary_dilation(cl20>0.03,iterations=2)&land&(cf>=0.5); gainm=binary_dilation(cg>0.03,iterations=2)&land
chg=np.full((ny,nx),np.nan); chg[land&(cf<0.5)]=0; chg[cf>=0.5]=1; chg[gainm]=3; chg[lossm]=2
axD.imshow(chg,extent=ext,origin="upper",cmap=ListedColormap([NF,STB,LS,GN]),vmin=-0.5,vmax=3.5); basemap(axD)
axD.set_title("(D) Forest change 2000–2020",fontsize=11,loc="left",fontweight="bold")
axE=fig.add_subplot(gs[1,1]); y=np.arange(len(bands))
axE.barh(y,decl,color="#e08214",edgecolor="#7a4a00",lw=.5,zorder=3)
for yi,d in zip(y,decl): axE.text(d+0.03,yi,f"{d:.2f}",va="center",ha="left",fontsize=8)
axE.set_yticks(y); axE.set_yticklabels(blab,fontsize=7.5); axE.set_ylabel("Elevation (m)",fontsize=9.5,labelpad=20)
axE.set_xlim(0,max(decl)*1.30); axE.set_xlabel("Forest decline 2000–2020 (% of band forest)",fontsize=8.6)
axE.xaxis.set_minor_locator(AutoMinorLocator()); axE.grid(axis="x",lw=0.4,color="0.85",zorder=0); axE.set_axisbelow(True)
axE.set_title("(E) Change by elevation",fontsize=11,loc="left",fontweight="bold"); axE.tick_params(labelsize=8)
axF=fig.add_subplot(gs[1,2]); yr=[2000,2010,2020]
axF.plot(yr,[a00,a10,a20],"-o",color=FG,lw=2.0,ms=6,zorder=3)
for xi,v,ha in zip(yr,[a00,a10,a20],["left","center","right"]):
    ox={"left":8,"center":0,"right":-8}[ha]; axF.annotate(f"{v:,.0f}",(xi,v),textcoords="offset points",xytext=(ox,9),fontsize=8.5,ha=ha,color=FG)
axF.set_ylabel("Forest area (km²)",fontsize=9.5,labelpad=22); axF.set_xlabel("Year",fontsize=9.5); axF.set_xticks(yr); axF.set_xlim(1997,2023)
axF.set_ylim(a20-500,a00+520); axF.yaxis.set_minor_locator(AutoMinorLocator()); axF.grid(axis="y",lw=0.4,color="0.85",zorder=0); axF.set_axisbelow(True)
axF.set_title("(F) Forest-area trend",fontsize=11,loc="left",fontweight="bold"); axF.tick_params(labelsize=8)
gsR=GridSpecFromSubplotSpec(2,1,subplot_spec=gs[:,3],height_ratios=[0.13,0.87],hspace=0.05)
axLeg=fig.add_subplot(gsR[0]); axLeg.axis("off")
axLeg.legend(handles=[Patch(fc=FG,label="Forest"),Patch(fc=NF,ec="0.5",lw=.3,label="Non-forest"),
             Patch(fc=STB,label="Stable forest"),Patch(fc=LS,label="Loss"),Patch(fc=GN,label="Gain")],
             loc="center",ncol=2,fontsize=9,frameon=False,handlelength=1.3,columnspacing=1.4,labelspacing=0.5)
axG=fig.add_subplot(gsR[1]); yb=np.arange(len(bands)); bl=[b.replace("\n","") for b in blab]
axG.barh(yb,[-x for x in loss_b],color=LS,edgecolor="#7a0000",lw=.4,zorder=3)
axG.barh(yb,gain_b,color=GN,edgecolor="#0d3b66",lw=.4,zorder=3)
axG.plot(net_b,yb,"-o",color="black",lw=1.3,ms=4,zorder=5)
for yi,nv in zip(yb,net_b): axG.text(nv,yi+0.30,f"{nv:+.0f}",fontsize=7.5,ha="center",va="bottom",color="black")
axG.axvline(0,color="0.3",lw=0.8,zorder=2)
axG.set_yticks(yb); axG.set_yticklabels(bl,fontsize=8); axG.set_ylabel("Elevation band (m)",fontsize=9,labelpad=2)
axG.set_xlabel("Forest-cover change 2000–2020 (km²)",fontsize=9)
axG.set_title("(G) Forest-cover change by elevation band",fontsize=10.5,loc="left",fontweight="bold")
axG.grid(axis="x",lw=0.4,color="0.85",zorder=0); axG.set_axisbelow(True); axG.tick_params(labelsize=8)
axG.legend(handles=[Patch(fc=LS,label="Loss"),Patch(fc=GN,label="Gain"),Line2D([0],[0],color="black",marker="o",ms=4,label="Net change")],
           loc="lower left",fontsize=8,frameon=True,handlelength=1.4)
# make (E) and (F) the same height as the maps (aspect-limited): match map (D) drawn box
fig.canvas.draw(); pd=axD.get_position()
for ax in (axE,axF):
    p=ax.get_position(); ax.set_position([p.x0,pd.y0,p.width,pd.height])
fig.suptitle("Two decades of forest dynamics — Eastern Black Sea Region (Türkiye), 2000–2020",fontsize=13.5,fontweight="bold",y=0.975)
fig.text(0.5,0.095,f"Forest area — 2000: {a00:,.0f} · 2010: {a10:,.0f} · 2020: {a20:,.0f} km².   Gross loss (2001–2020): {lo20:,.0f} km²; regrowth (2000–2012): {gain:,.0f} km²; net change: {net:,.0f} km² ({pnet:+.1f}%); canopy \u2265{TH}%.",ha="center",fontsize=9,color="0.12")
fig.text(0.5,0.058,"Forest extent is essentially stable; change is localized to valley sides and road corridors, with no broad decline. Loss/gain pixels in (D) are enlarged for visibility; true areas as stated above.",ha="center",fontsize=9,color="0.12")
fig.text(0.5,0.022,"Data: Hansen/UMD Global Forest Change v1.11 (30 m; tree cover \u2265%d%%, loss 2001–2020, gain 2000–2012). DEM: SRTM/GEBCO. Coast: Natural Earth. WGS84. Source: authors."%TH,ha="center",fontsize=8,style="italic",color="0.4")
fig.savefig("forest_dynamics.png",dpi=200,facecolor="white"); print("saved")
PY
