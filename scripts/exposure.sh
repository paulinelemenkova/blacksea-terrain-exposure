#!/usr/bin/env bash
# Combined figure - Conditioning-factor predisposition + population exposure, Eastern Black Sea.
# Compact layout (scientific-plotting): maps fill; predisposition legend + shared population
# colourbar below the maps; (d) at map height. (a) predisposition index; (b) population density
# (WorldPop 2020, 1 km); (c) population on high/very-high predisposition; (d) population by class.
set -e
REG=-R37.2/42.3/39.6/41.75
gmt grdcut @earth_relief_03s $REG -Gdem03.nc
gmt grdsample dem03.nc -I0.0025 -Gdem.nc
gmt grdconvert dem.nc dem.tif=gd:GTiff
curl -sL -o coast.geojson https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/geojson/ne_10m_coastline.geojson
python3 << 'PY'
import os,warnings; warnings.filterwarnings("ignore")
os.environ["GDAL_DISABLE_READDIR_ON_OPEN"]="EMPTY_DIR"
import numpy as np, rasterio, geopandas as gpd, requests, heapq
from rasterio.warp import reproject, Resampling
from rasterio.windows import from_bounds, bounds as wbounds
from rasterio.transform import from_bounds as tfb
from scipy.ndimage import binary_dilation
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import ListedColormap, LogNorm
from matplotlib.patches import Patch
from matplotlib.ticker import MultipleLocator, FuncFormatter, AutoMinorLocator
plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Nimbus Sans","Arial","Liberation Sans","DejaVu Sans"],"font.size":9})
W,E,S,N=37.2,42.3,39.6,41.75
dd=rasterio.open("dem.tif"); dem=dd.read(1).astype("float32"); dem[dem<-1e4]=np.nan
ny,nx=dem.shape; tr=dd.transform; crs=dd.crs; ext=[W,E,S,N]
latc=0.5*(S+N); dx=dd.res[0]*111320*np.cos(np.radians(latc)); dy=dd.res[1]*110540
land=dem>=0; sea=~land
gy=-np.gradient(dem,dy,axis=0); gx=np.gradient(dem,dx,axis=1); slope=np.degrees(np.arctan(np.hypot(gx,gy))); L=np.sqrt(dx*dy)
def sh(a,dr,dc):
    b=np.full_like(a,np.nan); sr=slice(max(dr,0),ny+min(dr,0)); srs=slice(max(-dr,0),ny+min(-dr,0)); sc=slice(max(dc,0),nx+min(dc,0)); scs=slice(max(-dc,0),nx+min(-dc,0)); b[srs,scs]=a[sr,sc]; return b
z2=sh(dem,1,0);z4=sh(dem,0,-1);z6=sh(dem,0,1);z8=sh(dem,-1,0);z1=sh(dem,1,-1);z3=sh(dem,1,1);z7=sh(dem,-1,-1);z9=sh(dem,-1,1)
D=((z4+z6)/2-dem)/L**2; Ec=((z2+z8)/2-dem)/L**2; Fh=(-z1+z3+z7-z9)/(4*L**2); G=(-z4+z6)/(2*L); H=(z2-z8)/(2*L)
den=np.where((G**2+H**2)==0,np.nan,G**2+H**2); plan=2*(D*H**2+Ec*G**2-Fh*G*H)/den
coast_sea=sea&binary_dilation(land); demf=np.where(sea,0.0,1e12).astype("float64"); vis=sea.copy(); heap=[]
for i,j in zip(*np.nonzero(coast_sea)): heapq.heappush(heap,(0.0,int(i),int(j)))
for i in range(ny):
    for j in (0,nx-1):
        if land[i,j] and not vis[i,j]: vis[i,j]=True; demf[i,j]=dem[i,j]; heapq.heappush(heap,(float(dem[i,j]),i,j))
for j in range(nx):
    for i in (0,ny-1):
        if land[i,j] and not vis[i,j]: vis[i,j]=True; demf[i,j]=dem[i,j]; heapq.heappush(heap,(float(dem[i,j]),i,j))
eps=1e-3; nbb=((-1,0),(1,0),(0,-1),(0,1),(-1,-1),(-1,1),(1,-1),(1,1))
while heap:
    e,i,j=heapq.heappop(heap)
    for di,dj in nbb:
        ni,nj=i+di,j+dj
        if 0<=ni<ny and 0<=nj<nx and not vis[ni,nj]:
            vis[ni,nj]=True; ne=dem[ni,nj]; ne=ne if ne>=e+eps else e+eps; demf[ni,nj]=ne; heapq.heappush(heap,(ne,ni,nj))
dw=np.where(sea,1e12,demf); offs=[(-1,-1),(-1,0),(-1,1),(0,-1),(0,1),(1,-1),(1,0),(1,1)]
dists=[np.hypot(dx,dy),dy,np.hypot(dx,dy),dx,dx,np.hypot(dx,dy),dy,np.hypot(dx,dy)]
recv=np.full(ny*nx,-1,np.int64); best=np.zeros((ny,nx)); rr,cc=np.mgrid[0:ny,0:nx]
for (dr,dc),ddd in zip(offs,dists):
    zn=sh(dw,dr,dc); drop=(dw-zn)/ddd; nr=rr+dr; ncc=cc+dc
    val=(nr>=0)&(nr<ny)&(ncc>=0)&(ncc<nx)&(drop>best)&np.isfinite(drop)&land
    recv.reshape(ny,nx)[val]=(nr*nx+ncc)[val]; best[val]=drop[val]
acc=np.ones(ny*nx); order=np.argsort(dw.ravel())[::-1]
for c in order:
    if recv[c]>=0: acc[recv[c]]+=acc[c]
acc=acc.reshape(ny,nx); acc[sea]=np.nan; As=acc*dx*dy/np.sqrt(dx*dy); tanb=np.tan(np.radians(np.maximum(slope,0.5)))
TWI=np.log(As/tanb); SPI=np.log(As*tanb+1.0)
def nz(x):
    v=x.copy().astype("float32"); v[sea]=np.nan; lo,hi=np.nanpercentile(v,2),np.nanpercentile(v,98); return np.clip((v-lo)/(hi-lo+1e-9),0,1)
comp=0.35*nz(slope)+0.25*nz(TWI)+0.20*nz(-plan)+0.20*nz(SPI); comp[sea]=np.nan
q=np.nanpercentile(comp,[20,40,60,80]); cls=np.digitize(comp,q).astype("float32"); cls[sea]=np.nan
pop1=None
for u in ["https://data.worldpop.org/GIS/Population/Global_2000_2020_1km/2020/TUR/tur_ppp_2020_1km_Aggregated.tif",
          "https://data.worldpop.org/GIS/Population/Global_2000_2020_1km_UNadj/2020/TUR/tur_ppp_2020_1km_Aggregated_UNadj.tif"]:
    try:
        r=requests.get(u,timeout=(15,180)); r.raise_for_status(); open("/tmp/pop.tif","wb").write(r.content)
        dp=rasterio.open("/tmp/pop.tif"); win=from_bounds(W,S,E,N,dp.transform)
        pop1=dp.read(1,window=win).astype("float32"); pop1[pop1<0]=0
        bb=wbounds(win,dp.transform); pop_wt=tfb(bb[0],bb[1],bb[2],bb[3],pop1.shape[1],pop1.shape[0]); break
    except Exception as ex: print("pop fail",str(ex)[:50])
clsi=np.where(sea,-1,cls).astype("float32")
cls_pop=np.full(pop1.shape,-1,"float32"); reproject(clsi,cls_pop,src_transform=tr,src_crs=crs,dst_transform=pop_wt,dst_crs=crs,resampling=Resampling.nearest,dst_nodata=-1)
popf=np.zeros((ny,nx),"float32"); reproject(pop1,popf,src_transform=pop_wt,src_crs=crs,dst_transform=tr,dst_crs=crs,resampling=Resampling.bilinear); popf[sea]=np.nan
total=float(pop1[cls_pop>=0].sum()); expc=[float(pop1[(cls_pop==k)].sum()) for k in range(5)]; share_hi=100*(expc[3]+expc[4])/total
coast=gpd.read_file("coast.geojson",bbox=(W,S,E,N)); AR=1/np.cos(np.radians(latc))
CIT=[("Giresun",38.39,40.91),("Trabzon",39.72,41.00),("Rize",40.52,41.02),("Artvin",41.82,41.18),("Ordu",37.88,40.98)]
cmap5=ListedColormap(["#1a9850","#a6d96a","#fee08b","#f46d43","#a50026"]); LAB=["Very low","Low","Moderate","High","Very high"]
fig=plt.figure(figsize=(16.5,3.9))
gs=GridSpec(1,4,figure=fig,width_ratios=[1,1,1,0.70],left=0.035,right=0.99,top=0.865,bottom=0.28,wspace=0.10)
def base(ax):
    coast.plot(ax=ax,color="black",lw=0.5,zorder=5); ax.set_xlim(W,E); ax.set_ylim(S,N); ax.set_aspect(AR)
    for xln in range(38,43): ax.axvline(xln,color="white",lw=0.4,alpha=0.6,zorder=4)
    for yln in range(40,42): ax.axhline(yln,color="white",lw=0.4,alpha=0.6,zorder=4)
    for n,lo,la in CIT: ax.plot(lo,la,"o",ms=3.2,mfc="black",mec="white",mew=0.5,zorder=6); ax.annotate(n,(lo,la),xytext=(0,4),textcoords="offset points",ha="center",fontsize=6.8,zorder=7)
    ax.xaxis.set_major_locator(MultipleLocator(1)); ax.yaxis.set_major_locator(MultipleLocator(1))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x,_:f"{x:g}°E")); ax.yaxis.set_major_formatter(FuncFormatter(lambda x,_:f"{x:g}°N")); ax.tick_params(labelsize=7.2)
axa=fig.add_subplot(gs[0,0]); axa.set_facecolor("#DCEAF4")
axa.imshow(cls,extent=ext,origin="upper",cmap=cmap5,vmin=-0.5,vmax=4.5); base(axa)
axa.set_title("(a) Predisposition index",fontsize=10,loc="left",fontweight="bold")
vmax=np.nanpercentile(popf[popf>0],99.5); NORM=LogNorm(vmin=1,vmax=vmax)
axb=fig.add_subplot(gs[0,1]); axb.set_facecolor("#DCEAF4")
imb=axb.imshow(np.where(popf>0.5,popf,np.nan),extent=ext,origin="upper",cmap="YlOrRd",norm=NORM); base(axb)
axb.set_title("(b) Population density (WorldPop 2020)",fontsize=10,loc="left",fontweight="bold")
axc=fig.add_subplot(gs[0,2]); axc.set_facecolor("#DCEAF4")
axc.imshow(np.where(land,0.9,np.nan),extent=ext,origin="upper",cmap=ListedColormap(["#eeeeee"]),vmin=0,vmax=1,zorder=0)
imc=axc.imshow(np.where((cls>=3)&(popf>0.5),popf,np.nan),extent=ext,origin="upper",cmap="YlOrRd",norm=NORM,zorder=1); base(axc)
axc.set_title("(c) Population on high / very-high predisposition",fontsize=9.6,loc="left",fontweight="bold")
axd=fig.add_subplot(gs[0,3]); yb=np.arange(5)
axd.barh(yb,[e/1000 for e in expc],height=0.66,color=cmap5.colors,edgecolor="0.3",lw=.4,zorder=3)
for yi,e in zip(yb,expc): axd.text(e/1000+total/1e5,yi,f"{e/1000:.0f}",va="center",ha="left",fontsize=7.5)
axd.set_yticks(yb); axd.set_yticklabels(LAB,fontsize=7.5)
axd.set_xlabel("Population (thousands)",fontsize=8.5); axd.set_title("(d) Population by class",fontsize=10,loc="left",fontweight="bold")
axd.xaxis.set_minor_locator(AutoMinorLocator()); axd.grid(axis="x",lw=0.4,color="0.85",zorder=0); axd.set_axisbelow(True); axd.tick_params(labelsize=7.5)
axd.text(0.98,0.97,f"≈ {total/1e6:.1f} M people\n{share_hi:.0f}% on High/Very-high\npredisposition terrain",transform=axd.transAxes,ha="right",va="top",fontsize=7.8,bbox=dict(fc="white",ec="0.6",lw=0.5,boxstyle="round,pad=0.3"))
fig.canvas.draw(); pa=axa.get_position(); pb=axb.get_position(); pc=axc.get_position(); pD=axd.get_position()
axd.set_position([pD.x0,pa.y0,pD.width,pa.height])
cax=fig.add_axes([pb.x0,pa.y0-0.11,pc.x1-pb.x0,0.032]); cb=fig.colorbar(imb,cax=cax,orientation="horizontal")
cb.set_label("Population density (persons km$^{-2}$)",fontsize=8); cb.ax.tick_params(labelsize=7)
fig.legend(handles=[Patch(fc=c,label=l) for c,l in zip(cmap5.colors,LAB)],loc="upper center",
           bbox_to_anchor=(pa.x0+pa.width/2,pa.y0-0.015),ncol=3,fontsize=7.5,frameon=True,
           title="Relative predisposition (a, c)",title_fontsize=8,handlelength=1.2,columnspacing=1.0,handletextpad=0.4)
fig.suptitle("Conditioning-factor predisposition and population exposure — Eastern Black Sea Region (Türkiye)",fontsize=12.5,fontweight="bold",y=0.955)
fig.text(0.5,0.03,"Heuristic predisposition (open-data overlay of slope, curvature, TWI, stream power; unvalidated, no inventory). Population: WorldPop 2020 (~1 km). Relief: SRTM/GEBCO. Coast: Natural Earth. WGS84.",ha="center",fontsize=7.5,style="italic",color="0.4")
fig.savefig("exposure.png",dpi=200,facecolor="white"); print("saved: pop=%.2fM share_hi=%.0f%%"%(total/1e6,share_hi))
PY
