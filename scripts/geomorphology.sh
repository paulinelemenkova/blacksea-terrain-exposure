#!/usr/bin/env bash
# Figure 4 - Geomorphology of the Eastern Black Sea Region (Türkiye)
# DEM from GMT earth_relief (SRTM/GEBCO), resampled ~180 m. Derivatives:
# slope, aspect (9-class), plan/profile curvature (Zevenbergen-Thorne),
# TRI (Riley), TWI (D8 flow accumulation), hypsometric curve, statistics.
# Needs: GMT 6.x, GDAL, python3 (numpy, scipy, matplotlib, rasterio, geopandas).
set -e
REG=-R37.2/42.3/39.6/41.75
gmt grdcut @earth_relief_03s $REG -Gdem03.nc          # or grdcut your GEBCO_2026.nc
gmt grdsample dem03.nc -I0.002 -Gdem.nc
gmt grdconvert dem.nc dem.tif=gd:GTiff
curl -sL -o coast.geojson https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/geojson/ne_10m_coastline.geojson

python3 << 'PY'
import warnings; warnings.filterwarnings("ignore")
import numpy as np, rasterio, geopandas as gpd, os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.patches import Patch
from matplotlib.ticker import MultipleLocator, FuncFormatter
W,E,S,N=37.2,42.3,39.6,41.75
ds=rasterio.open("dem.tif"); Z=ds.read(1).astype("float32"); Z[Z<-1e4]=np.nan
dem=np.where(Z>=0.0,Z,np.nan); ny,nx=dem.shape
dlon=(E-W)/nx; dlat=(N-S)/ny; latc=0.5*(S+N)
dx=dlon*111320*np.cos(np.radians(latc)); dy=dlat*110540; ext=[W,E,S,N]
gy=-np.gradient(dem,dy,axis=0); gx=np.gradient(dem,dx,axis=1)
slope=np.degrees(np.arctan(np.hypot(gx,gy)))
asp=np.degrees(np.arctan2(gy,-gx)); asp=np.where(asp<0,90-asp,np.where(asp>90,450-asp,90-asp))
flat=(np.hypot(gx,gy)<1e-6)|np.isnan(dem)
sect=((asp+22.5)//45).astype("int")%8
aspcls=np.where(flat,-1,sect).astype("float"); aspcls[np.isnan(dem)]=np.nan
L=np.sqrt(dx*dy)
def sh(a,dr,dc):
    b=np.full_like(a,np.nan)
    sr=slice(max(dr,0),ny+min(dr,0)); srs=slice(max(-dr,0),ny+min(-dr,0))
    sc=slice(max(dc,0),nx+min(dc,0)); scs=slice(max(-dc,0),nx+min(-dc,0))
    b[srs,scs]=a[sr,sc]; return b
z1=sh(dem,1,-1);z2=sh(dem,1,0);z3=sh(dem,1,1);z4=sh(dem,0,-1);z5=dem;z6=sh(dem,0,1)
z7=sh(dem,-1,-1);z8=sh(dem,-1,0);z9=sh(dem,-1,1)
D=((z4+z6)/2-z5)/L**2; Ec=((z2+z8)/2-z5)/L**2; F=(-z1+z3+z7-z9)/(4*L**2)
G=(-z4+z6)/(2*L); H=(z2-z8)/(2*L); den=np.where((G**2+H**2)==0,np.nan,G**2+H**2)
prof=-2*(D*G**2+Ec*H**2+F*G*H)/den; plan=2*(D*H**2+Ec*G**2-F*G*H)/den
nb=[sh(dem,dr,dc) for dr in(-1,0,1) for dc in(-1,0,1) if not(dr==0 and dc==0)]
tri=np.sqrt(np.nansum([(n-dem)**2 for n in nb],axis=0)); tri[np.isnan(dem)]=np.nan
zf=np.where(np.isnan(dem),1e9,dem)
offs=[(-1,-1),(-1,0),(-1,1),(0,-1),(0,1),(1,-1),(1,0),(1,1)]
dists=[np.hypot(dx,dy),dy,np.hypot(dx,dy),dx,dx,np.hypot(dx,dy),dy,np.hypot(dx,dy)]
recv=np.full(ny*nx,-1,dtype=np.int64); best=np.zeros((ny,nx),dtype="float32")
rr,cc=np.mgrid[0:ny,0:nx]
for (dr,dc),dd in zip(offs,dists):
    zn=sh(zf,-dr,-dc); drop=(zf-zn)/dd
    nr=rr+dr; ncc=cc+dc
    valid=(nr>=0)&(nr<ny)&(ncc>=0)&(ncc<nx)&(drop>best)&np.isfinite(drop)
    recv.reshape(ny,nx)[valid]=(nr*nx+ncc)[valid]; best[valid]=drop[valid]
acc=np.ones(ny*nx); order=np.argsort(zf.ravel())[::-1]
for c in order:
    rc=recv[c]
    if rc>=0: acc[rc]+=acc[c]
a_sca=acc.reshape(ny,nx)*dx*dy/np.sqrt(dx*dy)
twi=np.log(a_sca/np.tan(np.radians(np.maximum(slope,0.5)))); twi[np.isnan(dem)]=np.nan
v=dem[np.isfinite(dem)]
st=dict(Minimum=v.min(),Maximum=v.max(),Mean=v.mean(),Median=np.median(v),Std=v.std(),Relief=v.max()-v.min())
sv=np.sort(v)[::-1]; areafrac=np.arange(1,sv.size+1)/sv.size*100
coast=gpd.read_file("coast.geojson",bbox=(W,S,E,N))
CITIES=[("Ordu",37.88,40.98),("Giresun",38.39,40.91),("Trabzon",39.72,41.00),
        ("Rize",40.52,41.02),("Artvin",41.82,41.18),("Gumushane",39.48,40.46),("Bayburt",40.23,40.26)]
aspratio=1/np.cos(np.radians(latc)); plt.rcParams.update({"font.family":"DejaVu Sans","font.size":9})
fig=plt.figure(figsize=(16,11)); gs=GridSpec(3,4,hspace=0.32,wspace=0.30)
def mapfmt(ax,t):
    coast.plot(ax=ax,color="black",linewidth=0.5,zorder=5); ax.set_xlim(W,E); ax.set_ylim(S,N); ax.set_aspect(aspratio)
    ax.set_title(t,fontsize=10,loc="left"); ax.xaxis.set_major_locator(MultipleLocator(1)); ax.yaxis.set_major_locator(MultipleLocator(1))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x,_:f"{x:g}°E")); ax.yaxis.set_major_formatter(FuncFormatter(lambda x,_:f"{x:g}°N")); ax.tick_params(labelsize=7)
axA=fig.add_subplot(gs[0:2,0:2]); axA.set_facecolor("#DCEAF4")
im=axA.imshow(dem,extent=ext,origin="upper",cmap="terrain",vmin=0,vmax=3800); mapfmt(axA,"(A) Elevation (DEM ~180 m)")
for n,lo,la in CITIES:
    axA.plot(lo,la,"o",ms=4,mfc="black",mec="white",mew=0.6,zorder=6)
    axA.annotate(n,(lo,la),xytext=(0,4),textcoords="offset points",ha="center",fontsize=7,zorder=7)
fig.colorbar(im,ax=axA,fraction=0.03,pad=0.02).set_label("m",fontsize=8)
axB=fig.add_subplot(gs[0,2]); axB.set_facecolor("#DCEAF4")
cs=plt.get_cmap("YlOrRd"); imB=axB.imshow(slope,extent=ext,origin="upper",cmap=cs,norm=BoundaryNorm([0,5,15,25,35,45],cs.N,extend="max")); mapfmt(axB,"(B) Slope (°)")
fig.colorbar(imB,ax=axB,fraction=0.046,pad=0.02,extend="max")
axC=fig.add_subplot(gs[0,3]); axC.set_facecolor("#DCEAF4")
acolors=["#B0B0B0","#e41a1c","#ff7f00","#ffff33","#4daf4a","#00ced1","#377eb8","#984ea3","#f781bf"]; alab=["Flat","N","NE","E","SE","S","SW","W","NW"]
disp=np.where(aspcls==-1,0,aspcls+1); disp=np.where(np.isnan(aspcls),np.nan,disp)
imC=axC.imshow(disp,extent=ext,origin="upper",cmap=ListedColormap(acolors),vmin=0,vmax=9); mapfmt(axC,"(C) Aspect")
axC.legend(handles=[Patch(fc=acolors[i],ec="0.4",lw=0.2,label=alab[i]) for i in range(9)],fontsize=6,loc="center left",bbox_to_anchor=(1.01,0.5),frameon=True,handlelength=1)
axD=fig.add_subplot(gs[1,2]); axD.set_facecolor("#DCEAF4"); imD=axD.imshow(plan,extent=ext,origin="upper",cmap="RdBu",vmin=-0.02,vmax=0.02); mapfmt(axD,"(D) Plan curvature"); fig.colorbar(imD,ax=axD,fraction=0.046,pad=0.02,extend="both")
axE=fig.add_subplot(gs[1,3]); axE.set_facecolor("#DCEAF4"); imE=axE.imshow(prof,extent=ext,origin="upper",cmap="RdBu_r",vmin=-0.02,vmax=0.02); mapfmt(axE,"(E) Profile curvature"); fig.colorbar(imE,ax=axE,fraction=0.046,pad=0.02,extend="both")
axF=fig.add_subplot(gs[2,0]); axF.set_facecolor("#DCEAF4"); ct=plt.get_cmap("YlGnBu"); imF=axF.imshow(tri,extent=ext,origin="upper",cmap=ct,norm=BoundaryNorm([0,50,100,150,200],ct.N,extend="max")); mapfmt(axF,"(F) Terrain Ruggedness Index"); fig.colorbar(imF,ax=axF,fraction=0.046,pad=0.02,extend="max")
axG=fig.add_subplot(gs[2,1]); axG.set_facecolor("#DCEAF4"); imG=axG.imshow(twi,extent=ext,origin="upper",cmap="Blues",vmin=3,vmax=12); mapfmt(axG,"(G) Topographic Wetness Index"); fig.colorbar(imG,ax=axG,fraction=0.046,pad=0.02,extend="both")
axH=fig.add_subplot(gs[2,2]); axH.plot(areafrac,sv,color="#333",lw=1.5); axH.fill_between(areafrac,sv,color="#a6cee3",alpha=0.5)
axH.set_xlim(0,100); axH.set_ylim(0,st["Maximum"]*1.02); axH.set_xlabel("Area above elevation (%)"); axH.set_ylabel("Elevation (m)"); axH.set_title("(H) Hypsometric curve",fontsize=10,loc="left"); axH.grid(lw=0.4,color="0.85")
axI=fig.add_subplot(gs[2,3]); axI.axis("off"); axI.set_title("(I) Elevation statistics",fontsize=10,loc="left")
rows=[["Minimum (m)",f"{st['Minimum']:.0f}"],["Maximum (m)",f"{st['Maximum']:.0f}"],["Mean (m)",f"{st['Mean']:.0f}"],["Median (m)",f"{st['Median']:.0f}"],["Std. dev. (m)",f"{st['Std']:.0f}"],["Relief (m)",f"{st['Relief']:.0f}"]]
tb=axI.table(cellText=rows,colLabels=["Statistic","Value"],loc="center",cellLoc="left"); tb.auto_set_font_size(False); tb.set_fontsize(9); tb.scale(1,1.6)
fig.suptitle("Geomorphology of the Eastern Black Sea Region (Türkiye)",fontsize=15,fontweight="bold",y=0.985)
fig.text(0.5,0.02,"DEM: SRTM/GEBCO (GMT earth_relief), ~180 m. Coast: Natural Earth. WGS84.",ha="center",fontsize=8,style="italic",color="0.4")
fig.savefig("geomorphology.png",dpi=200,bbox_inches="tight",facecolor="white"); print("saved")
PY
