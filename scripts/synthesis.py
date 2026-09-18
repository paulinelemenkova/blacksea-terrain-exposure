#!/usr/bin/env python3
# Figure 8 - Conditioning-factor synthesis (heuristic predisposition index).
# Open-data weighted overlay (slope, plan-curvature convergence, TWI, SPI; +WorldClim
# rainfall where reachable). Black Sea masked (land=dem>=0). NOT a validated
# susceptibility model; no landslide inventory used. Needs GMT + python stack.
import os,warnings,subprocess; warnings.filterwarnings("ignore")
os.environ["GDAL_DISABLE_READDIR_ON_OPEN"]="EMPTY_DIR"
import numpy as np, rasterio, geopandas as gpd, requests, zipfile, io, heapq
from rasterio.warp import reproject, Resampling
from rasterio.windows import from_bounds, bounds as wbounds
from scipy.ndimage import binary_dilation
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
from matplotlib.ticker import MultipleLocator, FuncFormatter
plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Nimbus Sans","Arial","Liberation Sans","DejaVu Sans"],"font.size":9})
W,E,S,N=37.2,42.3,39.6,41.75
for c in ["gmt grdcut @earth_relief_03s -R37.2/42.3/39.6/41.75 -Gdem03.nc",
          "gmt grdsample dem03.nc -I0.0025 -Gdem.nc","gmt grdconvert dem.nc dem.tif=gd:GTiff"]:
    subprocess.run(c,shell=True,check=True)
if not os.path.exists("coast.geojson"):
    open("coast.geojson","wb").write(requests.get("https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/geojson/ne_10m_coastline.geojson",timeout=90).content)
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
rain=None
try:
    r=requests.get("https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_10m_prec.zip",timeout=(15,120)); r.raise_for_status()
    zipfile.ZipFile(io.BytesIO(r.content)).extractall("/tmp/p"); ann=None; wt=None
    for m in range(1,13):
        with rasterio.open(f"/tmp/p/wc2.1_10m_prec_{m:02d}.tif") as ds:
            win=from_bounds(W,S,E,N,ds.transform); a=ds.read(1,window=win,masked=True).astype("float32").filled(0)
            if wt is None: b=wbounds(win,ds.transform); wt=rasterio.transform.from_bounds(b[0],b[1],b[2],b[3],a.shape[1],a.shape[0])
            ann=a if ann is None else ann+a
    rain=np.zeros((ny,nx),"float32"); reproject(ann,rain,src_transform=wt,src_crs=crs,dst_transform=tr,dst_crs=crs,resampling=Resampling.bilinear)
except Exception as ex:
    print("rainfall unavailable:",str(ex)[:60])
def nz(x):
    v=x.copy().astype("float32"); v[sea]=np.nan; lo,hi=np.nanpercentile(v,2),np.nanpercentile(v,98); return np.clip((v-lo)/(hi-lo+1e-9),0,1)
f_slope=nz(slope); f_conv=nz(-plan); f_twi=nz(TWI); f_spi=nz(SPI)
if rain is not None:
    f_rain=nz(rain); comp=0.30*f_slope+0.20*f_twi+0.20*f_rain+0.15*f_conv+0.15*f_spi
    wtxt="Slope 0.30 \u00b7 TWI 0.20 \u00b7 Rainfall 0.20\nConvergence 0.15 \u00b7 Stream power 0.15"
else:
    comp=0.35*f_slope+0.25*f_twi+0.20*f_conv+0.20*f_spi
    wtxt="Slope 0.35 \u00b7 TWI 0.25 \u00b7 Convergence 0.20\nStream power 0.20  (rainfall added when\nWorldClim reachable)"
comp[sea]=np.nan; q=np.nanpercentile(comp,[20,40,60,80]); cls=np.digitize(comp,q).astype("float32"); cls[sea]=np.nan
coast=gpd.read_file("coast.geojson",bbox=(W,S,E,N)); AR=1/np.cos(np.radians(latc))
CIT=[("Giresun",38.39,40.91),("Trabzon",39.72,41.00),("Rize",40.52,41.02),("Artvin",41.82,41.18),("G\u00fcm\u00fc\u015fhane",39.48,40.46)]
fig=plt.figure(figsize=(15,8.0))
gs=GridSpec(5,3,figure=fig,width_ratios=[1,1,0.55],left=0.05,right=0.99,top=0.9,bottom=0.085,hspace=0.34,wspace=0.02)
axM=fig.add_subplot(gs[:,0:2]); axM.set_facecolor("#DCEAF4")
cmap5=ListedColormap(["#1a9850","#a6d96a","#fee08b","#f46d43","#a50026"])
axM.imshow(cls,extent=ext,origin="upper",cmap=cmap5,vmin=-0.5,vmax=4.5)
coast.plot(ax=axM,color="black",lw=0.6,zorder=5); axM.set_xlim(W,E); axM.set_ylim(S,N); axM.set_aspect(AR)
for n,lo,la in CIT:
    axM.plot(lo,la,"o",ms=4,mfc="black",mec="white",mew=0.6,zorder=6); axM.annotate(n,(lo,la),xytext=(0,5),textcoords="offset points",ha="center",fontsize=8,zorder=7)
axM.set_title("(a) Conditioning-factor synthesis \u2014 heuristic predisposition index",fontsize=12,loc="left",fontweight="bold")
axM.xaxis.set_major_locator(MultipleLocator(1)); axM.yaxis.set_major_locator(MultipleLocator(1))
axM.xaxis.set_major_formatter(FuncFormatter(lambda x,_:f"{x:g}\u00b0E")); axM.yaxis.set_major_formatter(FuncFormatter(lambda x,_:f"{x:g}\u00b0N")); axM.tick_params(labelsize=8)
for xln in range(38,43): axM.axvline(xln,color="white",lw=0.4,alpha=0.7,zorder=4)
for yln in range(40,42): axM.axhline(yln,color="white",lw=0.4,alpha=0.7,zorder=4)
axM.legend(handles=[Patch(fc=c,label=l) for c,l in zip(cmap5.colors,["Very low","Low","Moderate","High","Very high"])],
           title="Relative predisposition\n(quantile classes)",fontsize=8.5,title_fontsize=9,loc="lower left",frameon=True)
def rng(a,div=False):
    v=a.copy(); v[sea]=np.nan
    if div: lim=np.nanpercentile(np.abs(v),98); return v,-lim,lim
    return v,np.nanpercentile(v,2),np.nanpercentile(v,98)
smaps=[("(b) Slope",slope,"YlOrRd",False),("(c) Convergence (\u2212plan curv.)",-plan,"BrBG",True),
       ("(d) Wetness (TWI)",TWI,"jet_r",False),("(e) Stream power (SPI)",SPI,"rainbow",False)]
for k,(name,raw,cm,div) in enumerate(smaps):
    ax=fig.add_subplot(gs[k,2]); ax.set_facecolor("#DCEAF4"); v,vmn,vmx=rng(raw,div)
    ax.imshow(v,extent=ext,origin="upper",cmap=cm,vmin=vmn,vmax=vmx); coast.plot(ax=ax,color="black",lw=0.4,zorder=5)
    ax.set_xlim(W,E); ax.set_ylim(S,N); ax.set_aspect(AR)
    ax.xaxis.set_major_locator(MultipleLocator(2)); ax.yaxis.set_major_locator(MultipleLocator(1))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x,_:f"{x:g}\u00b0E")); ax.yaxis.set_major_formatter(FuncFormatter(lambda x,_:f"{x:g}\u00b0N"))
    ax.tick_params(labelsize=6,length=2,pad=1); ax.set_title(name,fontsize=8.6,loc="left")
    cb=fig.colorbar(im,ax=ax,fraction=0.05,pad=0.03); cb.ax.tick_params(labelsize=6)
axT=fig.add_subplot(gs[4,2]); axT.axis("off")
axT.text(0,1.0,"Weights (illustrative, AHP-style)",fontsize=9,fontweight="bold",va="top",transform=axT.transAxes)
axT.text(0,0.62,wtxt+"\nFactors min\u2013max normalised (2\u201398%).",fontsize=8,va="top",transform=axT.transAxes,linespacing=1.5)
fig.suptitle("Conditioning-factor synthesis for slope instability \u2014 Eastern Black Sea Region (T\u00fcrkiye)",fontsize=13.5,fontweight="bold",y=0.965)
fig.text(0.5,0.035,"Open-data heuristic overlay: SRTM/GEBCO DEM (slope, curvature, TWI, SPI)"+(" + WorldClim rainfall" if rain is not None else "")+". Coast: Natural Earth. WGS84. Heuristic; unvalidated; no landslide inventory.",ha="center",fontsize=8,style="italic",color="0.4")
fig.savefig("synthesis.png",dpi=200,bbox_inches="tight",facecolor="white"); print("saved (rain=%s)"%(rain is not None))
