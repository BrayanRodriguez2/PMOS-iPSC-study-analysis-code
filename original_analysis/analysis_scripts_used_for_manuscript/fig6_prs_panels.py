import numpy as np, pandas as pd, matplotlib, sys
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image
from scipy import stats
OUT, PRSX, OLD6 = sys.argv[1:4]
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":8,"axes.spines.top":False,"axes.spines.right":False,
 "axes.linewidth":0.6,"xtick.major.width":0.6,"ytick.major.width":0.6,"axes.edgecolor":"#555","axes.labelcolor":"#222",
 "xtick.color":"#444","ytick.color":"#444","axes.titlesize":9,"savefig.dpi":300})
C={"PMOS":"#eb6834","CTRL":"#2a78d6"}; M={"PMOS":"o","CTRL":"s"}; INK="#222"
cs=pd.read_excel(PRSX,sheet_name="Cohort_scores"); d=pd.read_excel(PRSX,sheet_name="iPSC_donors")
m=cs[cs.in_main_129==True]; assert len(m)==129
z={g:m.loc[m.group==g,"PRS_V1_z"].values for g in ("CTRL","PMOS")}
print("n",{g:len(v) for g,v in z.items()})
U=stats.mannwhitneyu(z["PMOS"],z["CTRL"]); auc=U.statistic/(len(z["PMOS"])*len(z["CTRL"]))
print("AUC",round(auc,3),"p",U.pvalue)
# ---- panel E ----
W,H=(2072-1450)/300,(1412-700)/300
fig,ax=plt.subplots(figsize=(W,H)); rng=np.random.default_rng(5)
for i,g in enumerate(("CTRL","PMOS")):
    v=z[g]
    vp=ax.violinplot(v,positions=[i],widths=.75,showextrema=False)
    for b in vp["bodies"]: b.set_facecolor(C[g]); b.set_alpha(.18); b.set_edgecolor("none")
    ax.scatter(i+rng.uniform(-.13,.13,len(v)),v,s=6,marker=M[g],color=C[g],edgecolor="white",lw=.25,zorder=3)
    ax.hlines(np.median(v),i-.26,i+.26,color=INK,lw=1.1,zorder=4)
ax.set_xticks([0,1]); ax.set_xticklabels([f"Control\n(n = {len(z['CTRL'])})",f"PMOS\n(n = {len(z['PMOS'])})"])
ax.set_xlim(-.6,1.6); ax.set_ylabel("PCOS PRS (z)")
ax.set_title("OR 1.81/SD (1.22–2.70)\nAUC 0.65 (0.55–0.74)",fontsize=6.5,color="#444")
ax.text(-0.36,1.03,"E",transform=ax.transAxes,fontsize=11,fontweight="bold",va="bottom")
fig.subplots_adjust(left=0.27,right=0.97,bottom=0.15,top=0.82); fig.savefig(f"{OUT}/_E.png",dpi=300); plt.close()
# ---- panel F ----
fig,ax=plt.subplots(figsize=(2072/300,2.25))
allz=m.PRS_V1_z.values
ax.hist(allz,bins=np.arange(-2.5,3.75,0.25),color="#d9dbde",edgecolor="white",lw=.4)
ymax=ax.get_ylim()[1]
ds=d.dropna(subset=["PRS_V1_z"]).sort_values("PRS_V1_z").reset_index(drop=True)
lev=[1.05,1.30,1.55,1.80]
for k,r in ds.iterrows():
    g=r.group; y=ymax*lev[k%4]
    ax.plot([r.PRS_V1_z]*2,[0,y],color=C[g],lw=.7,zorder=2)
    ax.scatter(r.PRS_V1_z,y,s=16,marker=M[g],color=C[g],edgecolor="white",lw=.4,zorder=3)
    ax.text(r.PRS_V1_z,y+ymax*.06,r.iPSC_line,ha="center",va="bottom",fontsize=5.6,color=C[g])
ax.set_ylim(0,ymax*2.1); ax.set_xlim(-2.6,3.7)
ax.set_xlabel("PCOS PRS (z within 129 women)"); ax.set_ylabel("Women")
mp=d.loc[d.group=="PMOS","PRS_V1_percentile"].mean(); mc=d.loc[d.group=="CTRL","PRS_V1_percentile"].mean()
ax.text(0.995,0.97,f"Mean percentile\nPMOS {mp:.0f} (n = 5)\ncontrol {mc:.0f} (n = 7)\n\nNot scored\n(contaminated array DNA):\nHEL351.3, HEL360.2",
        transform=ax.transAxes,ha="right",va="top",fontsize=6,color="#555")
ax.text(-0.075,0.98,"F",transform=ax.transAxes,fontsize=11,fontweight="bold",va="bottom")
fig.subplots_adjust(left=0.085,right=0.985,bottom=0.2,top=0.9); fig.savefig(f"{OUT}/_F.png",dpi=300); plt.close()
# ---- compose ----
old=Image.open(OLD6).convert("RGB"); E=Image.open(f"{OUT}/_E.png").convert("RGB"); F=Image.open(f"{OUT}/_F.png").convert("RGB")
top=old.copy(); top.paste((255,255,255),(1450,700,2072,1412)); top.paste(E.resize((622,712)),(1450,700))
canvas=Image.new("RGB",(2072,1412+F.height+20),"white"); canvas.paste(top,(0,0)); canvas.paste(F.resize((2072,F.height)),(0,1412+20))
canvas.save(f"{OUT}/Figure6.png",dpi=(300,300)); canvas.save(f"{OUT}/Figure6.tiff",dpi=(300,300),compression="tiff_lzw")
print("size",canvas.size)
