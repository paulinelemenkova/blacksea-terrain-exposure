#!/usr/bin/env python3
# =============================================================================
# Figure 2 - Eastern Black Sea Region (Türkiye): a matched pair of maps from the
# MTA/MRE 1:500,000 Geology Map of Türkiye (Orr & Associates, 2002).
#   (a) GEOLOGICAL  map - polygons coloured by chronostratigraphic AGE
#   (b) LITHOLOGICAL map - polygons coloured by LITHOLOGICAL ASSOCIATION
# Both: faults, mineral occurrences, province boundaries, coast, cities.
# Reads MapInfo TAB directly via GeoPandas/pyogrio. English legends.
# =============================================================================
import warnings; warnings.filterwarnings("ignore")
import numpy as np, geopandas as gpd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator, FuncFormatter

DATA  = "data/geology/"
COAST = "data/ne_coast.geojson"
ADM1  = "data/ne_adm1.geojson"
BBOX  = (37.2, 39.6, 42.3, 41.75)

# ---- AGE scheme (chronostratigraphic, young -> old) -------------------------
AGE = {
 "Holocene":"Quaternary","Pleistocene":"Quaternary","Quaternary":"Quaternary",
 "Pliocene":"Neogene","Miocene":"Neogene","Lower Miocene":"Neogene","Upper Miocene":"Neogene",
 "Oligocene-Miocene":"Neogene","Neogene":"Neogene","Tertiary":"Tertiary, undiff.",
 "Eocene":"Paleogene","Middle Eocene":"Paleogene","Oligocene":"Paleogene",
 "Upper Cretaceous":"Cretaceous","Lower Cretaceous":"Cretaceous","Cretaceous":"Cretaceous",
 "Jurassic-Cretaceous":"Jurassic-Cretaceous","Jurassic":"Jurassic","Mesozoic":"Mesozoic, undiff.",
 "Paleozoic":"Paleozoic","Permo-Carboniferous":"Paleozoic",
 "Unknown":"Age undetermined","Not applicable":"Age undetermined"}
AGE_ORDER = ["Quaternary","Neogene","Tertiary, undiff.","Paleogene","Cretaceous",
             "Jurassic-Cretaceous","Jurassic","Mesozoic, undiff.","Paleozoic",
             "Age undetermined","Water"]
AGE_COLOR = {"Quaternary":"#FCF7A6","Neogene":"#FFE619","Tertiary, undiff.":"#FCC98B",
             "Paleogene":"#FD9A52","Cretaceous":"#7FC64E","Jurassic-Cretaceous":"#8FD3A8",
             "Jurassic":"#34B2C9","Mesozoic, undiff.":"#A6DBC6","Paleozoic":"#A67A3B",
             "Age undetermined":"#D2D2D2","Water":"#CDE4F0"}

# ---- LITHOLOGICAL scheme ----------------------------------------------------
LIT = {
 "Unconsolidated and Semiconsolidated":"Unconsolidated & semiconsolidated",
 "Sedimentary":"Sedimentary rocks","Volcaniclastic":"Volcaniclastic rocks",
 "Volcanic":"Volcanic, undifferentiated","Volcanic-Hypabyssal":"Volcanic, undifferentiated",
 "Acid Volcanic-Hypabyssal":"Acid volcanic & hypabyssal",
 "Intermediate Volcanic-Hypabyssal":"Intermediate volcanic & hypabyssal",
 "Basic Volcanic-Hypabyssal":"Basic volcanic & hypabyssal",
 "Acid to intermediate Intrusive":"Granitoids & intrusive rocks",
 "Intermediate Intrusive":"Granitoids & intrusive rocks","Basic Intrusive":"Granitoids & intrusive rocks",
 "Ophiolite":"Ophiolite & ultrabasic rocks","Ultrabasic":"Ophiolite & ultrabasic rocks",
 "Metamorphic":"Metamorphic rocks","Water":"Water"}
LIT_ORDER = ["Unconsolidated & semiconsolidated","Sedimentary rocks","Volcaniclastic rocks",
             "Volcanic, undifferentiated","Acid volcanic & hypabyssal",
             "Intermediate volcanic & hypabyssal","Basic volcanic & hypabyssal",
             "Granitoids & intrusive rocks","Ophiolite & ultrabasic rocks",
             "Metamorphic rocks","Water"]
LIT_COLOR = {"Unconsolidated & semiconsolidated":"#FDF3C0","Sedimentary rocks":"#D9C08C",
             "Volcaniclastic rocks":"#B8A45C","Volcanic, undifferentiated":"#E34A33",
             "Acid volcanic & hypabyssal":"#F4A6C0","Intermediate volcanic & hypabyssal":"#C2593F",
             "Basic volcanic & hypabyssal":"#7A5195","Granitoids & intrusive rocks":"#C5197D",
             "Ophiolite & ultrabasic rocks":"#2CA25F","Metamorphic rocks":"#9E7BB5","Water":"#CDE4F0"}

def box():
    return gpd.GeoSeries.from_wkt(
        [f"POLYGON(({BBOX[0]} {BBOX[1]},{BBOX[2]} {BBOX[1]},{BBOX[2]} {BBOX[3]},"
         f"{BBOX[0]} {BBOX[3]},{BBOX[0]} {BBOX[1]}))"], crs=4326).iloc[0]

# ---- load once --------------------------------------------------------------
lith = gpd.read_file(DATA+"Turkey 500k Lithology V1.TAB", bbox=BBOX)
lith["geometry"] = lith.geometry.buffer(0); lith = lith.clip(box())
faults = gpd.read_file(DATA+"Turkey 500k Faults V1.TAB", bbox=BBOX).clip(box())
minocc = gpd.read_file(DATA+"Turkey 500k MinOcc V1.TAB", bbox=BBOX)
coast  = gpd.read_file(COAST, bbox=BBOX).clip(box())
prov   = gpd.read_file(ADM1, bbox=BBOX)
prov   = prov[prov["adm0_a3"].eq("TUR")].clip(box())
CITIES=[("Ordu",37.88,40.98),("Giresun",38.39,40.91),("Trabzon",39.72,41.00),
        ("Rize",40.52,41.02),("Artvin",41.82,41.18),("Gumushane",39.48,40.46),
        ("Bayburt",40.23,40.26)]
ASPECT = 1.0/np.cos(np.radians(0.5*(BBOX[1]+BBOX[3])))
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":11})

def make_map(mode, out, title):
    if mode=="age":
        cls = lith["Age"].map(lambda v: AGE.get(v,"Age undetermined"))
        cls = cls.where(lith["Lith_Association"]!="Water","Water")
        order,color,ltitle = AGE_ORDER,AGE_COLOR,"Age"
    else:
        cls = lith["Lith_Association"].map(lambda v: LIT.get(v,"Water"))
        order,color,ltitle = LIT_ORDER,LIT_COLOR,"Lithological association"
    cols = cls.map(color)

    fig, ax = plt.subplots(figsize=(15,9)); ax.set_facecolor("#DCEAF4")
    lith.plot(ax=ax, color=cols, edgecolor="0.4", linewidth=0.12)
    prov.boundary.plot(ax=ax, color="0.15", linewidth=0.7, linestyle=(0,(6,2,1,2)), alpha=0.7)
    faults[faults["Type"]=="Probable"].plot(ax=ax,color="#7a0000",linewidth=0.6,linestyle=(0,(4,2)))
    faults[faults["Type"]=="Mapped"].plot(ax=ax,color="#7a0000",linewidth=0.9)
    coast.plot(ax=ax, color="black", linewidth=0.8)
    minocc.plot(ax=ax, marker="D", markersize=26, facecolor="#fffbcc", edgecolor="black",
                linewidth=0.6, zorder=6)
    for n,lo,la in CITIES:
        ax.plot(lo,la,"o",ms=6,mfc="black",mec="white",mew=0.8,zorder=7)
        ax.annotate(n,(lo,la),xytext=(0,-13),textcoords="offset points",ha="center",
                    fontsize=9.5,fontweight="bold",zorder=8,
                    bbox=dict(boxstyle="round,pad=0.12",fc="white",ec="none",alpha=0.72))
    ax.text(40.2,41.58,"BLACK SEA",style="italic",fontsize=15,color="#20507a",ha="center",va="center")
    ax.set_xlim(BBOX[0],BBOX[2]); ax.set_ylim(BBOX[1],BBOX[3]); ax.set_aspect(ASPECT)
    ax.xaxis.set_major_locator(MultipleLocator(1)); ax.yaxis.set_major_locator(MultipleLocator(1))
    ax.xaxis.set_minor_locator(MultipleLocator(0.5)); ax.yaxis.set_minor_locator(MultipleLocator(0.5))
    ax.grid(which="major", color="white", linewidth=0.5, alpha=0.7)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v,_:f"{v:g}°E"))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v,_:f"{v:g}°N"))
    ax.tick_params(top=True,right=True,labelright=True,direction="out")
    for s in ax.spines.values(): s.set_linewidth(1.2)
    dlon=100/(111.32*np.cos(np.radians(40.0))); x0,y0=BBOX[0]+0.25,BBOX[1]+0.14
    ax.plot([x0,x0+dlon],[y0,y0],color="black",lw=3,solid_capstyle="butt")
    ax.text(x0+dlon/2,y0+0.075,"100 km",ha="center",va="bottom",fontsize=9,
            bbox=dict(boxstyle="round,pad=0.15",fc="white",ec="none",alpha=0.75))
    ax.annotate("N",xy=(37.55,41.62),xytext=(37.55,41.40),ha="center",fontsize=13,
                fontweight="bold",arrowprops=dict(arrowstyle="-|>",color="black",lw=1.6))
    ax.set_title(title, fontsize=13, fontweight="bold", pad=8)
    ax.text(0.5,-0.085,"Lithology & faults: MTA/MRE 1:500,000 Geology Map of Türkiye "
            "(Orr & Associates, 2002). Coast & provinces: Natural Earth. WGS84.",
            transform=ax.transAxes, ha="center", fontsize=8, style="italic", color="0.35")
    present=set(cls)
    handles=[Patch(fc=color[c],ec="0.4",lw=0.3,label=c) for c in order if c in present]
    handles+=[Line2D([0],[0],color="#7a0000",lw=1.2,label="Fault (mapped)"),
              Line2D([0],[0],color="#7a0000",lw=1.0,ls=(0,(4,2)),label="Fault (probable)"),
              Line2D([0],[0],color="0.15",lw=0.9,ls=(0,(6,2,1,2)),label="Province boundary"),
              Line2D([0],[0],marker="D",color="none",mfc="#fffbcc",mec="black",mew=0.6,
                     ms=8,label="Mineral occurrence")]
    leg=ax.legend(handles=handles,loc="center left",bbox_to_anchor=(1.005,0.5),frameon=True,
                  fontsize=9,title=ltitle,title_fontsize=10,borderpad=0.8,labelspacing=0.5)
    leg.get_frame().set_edgecolor("0.3")
    plt.savefig(out,dpi=300,bbox_inches="tight",facecolor="white"); plt.close()
    print("saved",out)

make_map("age","figures/geological_setting.png",
         "Geological setting of the Eastern Black Sea Region (Türkiye) — units by age")
make_map("lithology","figures/lithological_setting.png",
         "Lithological setting of the Eastern Black Sea Region (Türkiye) — rock associations")
