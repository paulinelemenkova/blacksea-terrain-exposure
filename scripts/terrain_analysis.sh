#!/usr/bin/env bash
# Figure 7 - Terrain analysis, Eastern Black Sea Region (Türkiye)  [rev.2]
# NEW layers only (no overlap with fig04): SPI, STI, TPI, Weiss landform classes,
# local relief, drainage network (blue), NW-SE profile with hypsometric fill.
# Priority-flood depression filling; profile sampled on real elevation, cropped to land.
set -e
REG=-R37.2/42.3/39.6/41.75
gmt grdcut @earth_relief_03s $REG -Gdem03.nc
gmt grdsample dem03.nc -I0.003 -Gdem.nc
gmt grdconvert dem.nc dem.tif=gd:GTiff
curl -sL -o coast.geojson https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/geojson/ne_10m_coastline.geojson
python3 << 'PY'
import warnings; warnings.filterwarnings("ignore")
import numpy as np, rasterio, geopandas as gpd, heapq
from scipy.ndimage import uniform_filter, maximum_filter, minimum_filter, map_coordinates, binary_dilation
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import BoundaryNorm, ListedColormap, LightSource, LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.patches import Patch, Polygon
from matplotlib.ticker import MultipleLocator, FuncFormatter
plt.rcParams.update({"font.family":"sans-serif",
  "font.sans-serif":["Nimbus Sans","Arial","Liberation Sans","Arimo","Helvetica","DejaVu Sans"],"font.size":9})
W,E,S,N=37.2,42.3,39.6,41.75
ds=rasterio.open("dem.tif"); Z=ds.read(1).astype("float32"); Z[Z<-1e4]=np.nan
dem=np.where(Z>=0.0,Z,np.nan); ny,nx=dem.shape
dlon=(E-W)/nx; dlat=(N-S)/ny; latc=0.5*(S+N)
dx=dlon*111320*np.cos(np.radians(latc)); dy=dlat*110540; ext=[W,E,S,N]; cellkm2=dx*dy/1e6
gy=-np.gradient(dem,dy,axis=0); gx=np.gradient(dem,dx,axis=1)
slope=np.arctan(np.hypot(gx,gy)); slp_deg=np.degrees(slope)
sea=np.isnan(dem); land=~sea; coast_sea=sea&binary_dilation(land)
demf=np.where(sea,0.0,1e12).astype("float64"); visited=sea.copy(); heap=[]
for i,j in zip(*np.nonzero(coast_sea)): heapq.heappush(heap,(0.0,int(i),int(j)))
for i in range(ny):
    for j in (0,nx-1):
        if land[i,j] and not visited[i,j]: visited[i,j]=True; demf[i,j]=dem[i,j]; heapq.heappush(heap,(float(dem[i,j]),i,j))
for j in range(nx):
    for i in (0,ny-1):
        if land[i,j] and not visited[i,j]: visited[i,j]=True; demf[i,j]=dem[i,j]; heapq.heappush(heap,(float(dem[i,j]),i,j))
eps=1e-3; nbr=((-1,0),(1,0),(0,-1),(0,1),(-1,-1),(-1,1),(1,-1),(1,1))
while heap:
    e,i,j=heapq.heappop(heap)
    for di,dj in nbr:
        ni=i+di; nj=j+dj
        if 0<=ni<ny and 0<=nj<nx and not visited[ni,nj]:
            visited[ni,nj]=True; ne=dem[ni,nj]
            if ne<e+eps: ne=e+eps
            demf[ni,nj]=ne; heapq.heappush(heap,(ne,ni,nj))
def sh(a,dr,dc):
    b=np.full_like(a,np.nan)
    sr=slice(max(dr,0),ny+min(dr,0)); srs=slice(max(-dr,0),ny+min(-dr,0))
    sc_=slice(max(dc,0),nx+min(dc,0)); scs=slice(max(-dc,0),nx+min(-dc,0)); b[srs,scs]=a[sr,sc_]; return b
demfw=np.where(sea,1e12,demf)
offs=[(-1,-1),(-1,0),(-1,1),(0,-1),(0,1),(1,-1),(1,0),(1,1)]
dists=[np.hypot(dx,dy),dy,np.hypot(dx,dy),dx,dx,np.hypot(dx,dy),dy,np.hypot(dx,dy)]
recv=np.full(ny*nx,-1,dtype=np.int64); best=np.zeros((ny,nx),dtype="float64"); rr,cc=np.mgrid[0:ny,0:nx]
for (dr,dc),dd in zip(offs,dists):
    zn=sh(demfw,dr,dc); drop=(demfw-zn)/dd; nr=rr+dr; ncc=cc+dc
    valid=(nr>=0)&(nr<ny)&(ncc>=0)&(ncc<nx)&(drop>best)&np.isfinite(drop)&land
    recv.reshape(ny,nx)[valid]=(nr*nx+ncc)[valid]; best[valid]=drop[valid]
acc=np.ones(ny*nx); order=np.argsort(demfw.ravel())[::-1]
for c in order:
    rc=recv[c]
    if rc>=0: acc[rc]+=acc[c]
acc=acc.reshape(ny,nx); acc[sea]=np.nan
As=acc*dx*dy/np.sqrt(dx*dy)
tanb=np.tan(np.maximum(slope,np.radians(0.5))); sinb=np.sin(np.maximum(slope,np.radians(0.5)))
SPI=np.log(As*tanb+1.0); SPI[sea]=np.nan
STI=(As/22.13)**0.6*(sinb/0.0896)**1.3; STI[sea]=np.nan
k=7; dm=np.where(sea,np.nanmean(dem),dem)
TPI=dem-uniform_filter(dm,size=k); TPI[sea]=np.nan
LR=(maximum_filter(dm,size=k)-minimum_filter(dm,size=k)); LR[sea]=np.nan
sd=np.nanstd(TPI); lf=np.full((ny,nx),np.nan)
lf=np.where(TPI>sd,6,lf); lf=np.where((TPI>0.5*sd)&(TPI<=sd),5,lf)
lf=np.where((TPI>=-0.5*sd)&(TPI<=0.5*sd)&(slp_deg>5),4,lf); lf=np.where((TPI>=-0.5*sd)&(TPI<=0.5*sd)&(slp_deg<=5),3,lf)
lf=np.where((TPI>=-sd)&(TPI<-0.5*sd),2,lf); lf=np.where(TPI<-sd,1,lf); lf[sea]=np.nan
strm=(acc*cellkm2)>8.0; strm=strm&land; strm=binary_dilation(strm)
ls=LightSource(azdeg=315,altdeg=45); hs=ls.hillshade(dm,vert_exag=1,dx=dx,dy=dy); hs[sea]=np.nan
lon0,lat0,lon1,lat1=39.3,41.2,41.7,39.9; npv=600
lons=np.linspace(lon0,lon1,npv); lats=np.linspace(lat0,lat1,npv)
demz=np.where(sea,0.0,dem)
prof=map_coordinates(demz,[(N-lats)/dlat,(lons-W)/dlon],order=1)
lline=map_coordinates(land.astype(float),[(N-lats)/dlat,(lons-W)/dlon],order=1)>0.5
R=6371.0; seg=[R*np.hypot(np.radians((lat1-lat0)/(npv-1)),np.radians((lon1-lon0)/(npv-1))*np.cos(np.radians(lat0+(lat1-lat0)*(i-1)/(npv-1)))) for i in range(1,npv)]
dfull=np.concatenate([[0],np.cumsum(seg)]); idx=np.where(lline)[0]; a0,a1=idx[0],idx[-1]
prof=prof[a0:a1+1]; dist=dfull[a0:a1+1]-dfull[a0]
nwll=f"{lons[a0]:.3f}°E, {lats[a0]:.3f}°N"; sell=f"{lons[a1]:.3f}°E, {lats[a1]:.3f}°N"
hyps=LinearSegmentedColormap.from_list("hyps",[(p/3800.0,c) for p,c in
   [(0,"#1a9850"),(500,"#66bd63"),(1000,"#a6d96a"),(1500,"#fee08b"),(2000,"#fdae61"),(2500,"#b2823f"),(3000,"#cdbfb0"),(3800,"#ffffff")]],N=256)
coast=gpd.read_file("coast.geojson",bbox=(W,S,E,N)); AR=1/np.cos(np.radians(latc))
fig=plt.figure(figsize=(16,8.8))
gs=GridSpec(3,3,figure=fig,left=0.045,right=0.985,top=0.94,bottom=0.10,hspace=0.20,wspace=0.15,height_ratios=[1,1,0.55])
def mapfmt(ax,t):
    coast.plot(ax=ax,color="black",linewidth=0.5,zorder=5); ax.set_xlim(W,E); ax.set_ylim(S,N); ax.set_aspect(AR)
    ax.set_title(t,fontsize=10,loc="left"); ax.xaxis.set_major_locator(MultipleLocator(1)); ax.yaxis.set_major_locator(MultipleLocator(1))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x,_:f"{x:g}°E")); ax.yaxis.set_major_formatter(FuncFormatter(lambda x,_:f"{x:g}°N")); ax.tick_params(labelsize=7)
def cb(im,ax,label,**kw):
    c=fig.colorbar(im,ax=ax,fraction=0.046,pad=0.02,**kw); c.set_label(label,fontsize=7.5); c.ax.tick_params(labelsize=6.5); return c
def rob(a,lo=2,hi=98): return np.nanpercentile(a,lo),np.nanpercentile(a,hi)
axa=fig.add_subplot(gs[0,0]); axa.set_facecolor("#DCEAF4"); v0,v1=rob(SPI)
im=axa.imshow(SPI,extent=ext,origin="upper",cmap="viridis",vmin=v0,vmax=v1); mapfmt(axa,"(a) Stream Power Index (ln)"); cb(im,axa,"ln(SPI)",extend="both")
axb=fig.add_subplot(gs[0,1]); axb.set_facecolor("#DCEAF4")
vs=STI[np.isfinite(STI)]; b=np.unique(np.percentile(vs,np.linspace(0,100,257)))
cmapb=plt.get_cmap("turbo",len(b)-1); normb=BoundaryNorm(b,len(b)-1)
im=axb.imshow(STI,extent=ext,origin="upper",cmap=cmapb,norm=normb); mapfmt(axb,"(b) Sediment Transport Index")
tks=np.percentile(vs,[1,25,50,75,99]); cbb=fig.colorbar(im,ax=axb,fraction=0.046,pad=0.02,ticks=tks)
cbb.ax.set_yticklabels([f"{t:.0f}" for t in tks],fontsize=6.5); cbb.set_label("STI",fontsize=7.5)
axc=fig.add_subplot(gs[0,2]); axc.set_facecolor("#DCEAF4"); lim=np.nanpercentile(np.abs(TPI),98)
im=axc.imshow(TPI,extent=ext,origin="upper",cmap="RdBu_r",vmin=-lim,vmax=lim); mapfmt(axc,"(c) Topographic Position Index"); cb(im,axc,"TPI (m)",extend="both")
axd=fig.add_subplot(gs[1,0]); axd.set_facecolor("#DCEAF4"); v0,v1=rob(LR)
im=axd.imshow(LR,extent=ext,origin="upper",cmap="YlGnBu",vmin=0,vmax=v1); mapfmt(axd,"(d) Local relief (1 km window)"); cb(im,axd,"m",extend="max")
axd.plot(lons[a0:a1+1],lats[a0:a1+1],color="red",lw=2.4,solid_capstyle="round",zorder=8)
axd.plot([lons[a0],lons[a1]],[lats[a0],lats[a1]],"o",ms=3.5,color="red",mec="white",mew=0.6,zorder=9)
axd.annotate("NW",(lons[a0],lats[a0]),xytext=(4,3),textcoords="offset points",fontsize=7.5,fontweight="bold",color="red",zorder=9)
axd.annotate("SE",(lons[a1],lats[a1]),xytext=(-16,-2),textcoords="offset points",fontsize=7.5,fontweight="bold",color="red",zorder=9)
axe=fig.add_subplot(gs[1,1]); axe.set_facecolor("#DCEAF4")
lcol=["#2166ac","#67a9cf","#d9f0d3","#fddbc7","#ef8a62","#b2182b"]; llab=["Valley","Lower slope","Flat","Mid slope","Upper slope","Ridge"]
im=axe.imshow(lf,extent=ext,origin="upper",cmap=ListedColormap(lcol),norm=BoundaryNorm(np.arange(0.5,7.5,1),6)); mapfmt(axe,"(e) Landform classes (Weiss)")
axe.legend(handles=[Patch(fc=lcol[i],ec="0.4",lw=0.2,label=llab[i]) for i in range(6)],fontsize=6,loc="center left",bbox_to_anchor=(1.01,0.5),frameon=True,handlelength=1)
axf=fig.add_subplot(gs[1,2]); axf.set_facecolor("#DCEAF4")
axf.imshow(hs,extent=ext,origin="upper",cmap="gray",vmin=0.15,vmax=1,zorder=1)
axf.imshow(np.where(strm,1.0,np.nan),extent=ext,origin="upper",cmap=ListedColormap(["#1f6fe0"]),zorder=2)
mapfmt(axf,"(f) Drainage network (> 8 km$^2$)")
axg=fig.add_subplot(gs[2,:]); ymax=float(np.nanmax(prof))*1.05
img=axg.imshow(np.linspace(0,3800,300).reshape(-1,1),aspect="auto",cmap=hyps,vmin=0,vmax=3800,extent=[0,dist[-1],0,ymax],origin="lower",zorder=1)
clip=Polygon(np.column_stack([np.r_[0,dist,dist[-1]],np.r_[0,prof,0]]),closed=True,fc="none",ec="none",transform=axg.transData)
axg.add_patch(clip); img.set_clip_path(clip); axg.plot(dist,prof,color="#333",lw=0.8,zorder=3)
axg.set_xlim(0,dist[-1]); axg.set_ylim(0,ymax); axg.grid(lw=0.4,color="0.85",zorder=0)
axg.set_xticks(list(range(0,int(dist[-1]),50))+[int(round(dist[-1]))])
axg.set_xlabel("Distance along NW–SE transect (km)",fontsize=9); axg.set_ylabel("Elevation (m)",fontsize=9)
axg.set_title("(g) Topographic profile (NW → SE)",fontsize=10,loc="left")
axg.text(0.01,0.90,"NW (coast)",transform=axg.transAxes,fontsize=8,color="#222")
axg.text(0.99,0.90,"SE (interior)",transform=axg.transAxes,fontsize=8,color="#222",ha="right")
axg.text(0.008,0.05,"NW  "+nwll,transform=axg.transAxes,fontsize=7.6,color="#c81e1e",ha="left",va="bottom")
axg.text(0.992,0.05,"SE  "+sell,transform=axg.transAxes,fontsize=7.6,color="#c81e1e",ha="right",va="bottom")
sm=ScalarMappable(cmap=hyps,norm=Normalize(0,3800)); sm.set_array([])
cg=fig.colorbar(sm,ax=axg,fraction=0.012,pad=0.01); cg.set_label("Elevation (m)",fontsize=8); cg.ax.tick_params(labelsize=6.5)
fig.suptitle("Terrain analysis of the Eastern Black Sea Region (Türkiye)",fontsize=14,fontweight="bold",y=0.985)
fig.text(0.5,0.048,"Transect endpoints of the plotted land profile (red line in panel d):  NW / left (coast) "+nwll+";   SE / right (interior) "+sell+".",ha="center",fontsize=8.6,color="0.18")
fig.text(0.5,0.018,"Hydro-geomorphometric indices from SRTM/GEBCO DEM (~330 m; depression-filled flow routing). Coast: Natural Earth. WGS84.",ha="center",fontsize=8,style="italic",color="0.45")
fig.savefig("terrain_analysis.png",dpi=200,facecolor="white"); print("saved")
PY
