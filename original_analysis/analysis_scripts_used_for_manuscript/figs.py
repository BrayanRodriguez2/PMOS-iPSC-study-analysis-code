import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from scipy import stats
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":8,"axes.spines.top":False,"axes.spines.right":False,
  "axes.linewidth":0.6,"xtick.major.width":0.6,"ytick.major.width":0.6,"axes.edgecolor":"#555","axes.labelcolor":"#222",
  "xtick.color":"#444","ytick.color":"#444","axes.titlesize":9,"axes.titleweight":"bold","savefig.dpi":300})
O="/home/claude/ms"; F=O+"/figs"
C={"PMOS":"#eb6834","control":"#2a78d6"}; M={"PMOS":"o","control":"s"}; INK="#222"; MUTED="#888"
key=pd.read_csv("/home/claude/pheno/sample_key_14_ipsc_donors.csv",dtype={"array_iid":str}).set_index("hel_line")
lg=pd.read_csv(f"{O}/tables/log2_expressed.tsv.gz",sep="\t",index_col=0)
cnt=pd.read_csv("/mnt/user-data/attachments/522de043-d88c-536e-9f44-66e7182b8053/336ebae1-Normalized_counts.tabular",sep="\t",index_col=0)
cnt.columns=["HEL"+c.split("_")[0][3:]+"."+c.split("_")[1] for c in cnt.columns]
key=key.loc[lg.columns]
def lab(ax,s): ax.text(-0.18,1.06,s,transform=ax.transAxes,fontsize=11,fontweight="bold",va="bottom")
def strip(ax,col,ylabel,log=False):
    rng=np.random.default_rng(3)
    for i,g in enumerate(["control","PMOS"]):
        v=key.loc[key.group==g,col].values
        x=i+rng.uniform(-.12,.12,len(v))
        ax.scatter(x,v,s=26,marker=M[g],color=C[g],edgecolor="white",linewidth=.6,zorder=3)
        ax.hlines(np.median(v),i-.25,i+.25,color=INK,lw=1.2,zorder=4)
    p=stats.mannwhitneyu(key.loc[key.group=="PMOS",col],key.loc[key.group=="control",col]).pvalue
    ax.set_xticks([0,1]); ax.set_xticklabels(["Control","PMOS"]); ax.set_xlim(-.6,1.6); ax.set_ylabel(ylabel)
    if log: ax.set_yscale("log")
    ax.set_title(f"p = {p:.3f}" if p>=0.001 else f"p = {p:.1e}",fontsize=7.5,fontweight="normal",color="#444")
 
# ---------- Fig 1 design ----------
fig,ax=plt.subplots(figsize=(7.2,3.9)); ax.axis("off"); ax.set_xlim(0,120); ax.set_ylim(0,54)
def box(x,y,w,h,t,fc="#f3f4f6",ec="#9aa0a6"):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.3,rounding_size=1.2",fc=fc,ec=ec,lw=.8))
    ax.text(x+w/2,y+h/2,t,ha="center",va="center",fontsize=6.6,color=INK,linespacing=1.35)
def arr(x1,y1,x2,y2): ax.annotate("",xy=(x2,y2),xytext=(x1,y1),arrowprops=dict(arrowstyle="-|>",color="#666",lw=.9,shrinkA=0,shrinkB=0))
ax.text(14,52,"Donors",ha="center",fontsize=7.5,fontweight="bold"); ax.text(50,52,"Cells and data",ha="center",fontsize=7.5,fontweight="bold"); ax.text(96,52,"Analyses",ha="center",fontsize=7.5,fontweight="bold")
box(1,36,27,12,"Oulu PMOS/control cohort\n140 women genotyped\n14 selected as iPSC donors")
box(1,19,27,12,"Clinical phenotyping\nAMH, testosterone, FAI,\nFG score, BMI, HOMA-IR")
box(1,2,27,12,"Illumina GSA (ESTchip-2)\n777,016 markers\nGRCh38")
box(36,36,28,12,"Endometrial stromal\nfibroblasts (eSF)",fc="#fde9df",ec=C["PMOS"])
box(36,19,28,12,"CRISPRa reprogramming\n7 PMOS + 7 control iPSC\n(Khatun et al. 2026)",fc="#fde9df",ec=C["PMOS"])
box(36,2,28,12,"Bulk mRNA-seq,\nundifferentiated iPSC\n16,489 expressed genes",fc="#e3eefb",ec=C["control"])
box(72,36,47,12,"1  Transcriptome vs diagnosis, hormones,\n    BMI and age (exact / permutation tests)")
box(72,19,47,12,"2  XIST erosion as a confounder\n3  Expression of PCOS GWAS-locus genes")
box(72,2,47,12,"4  Targeted cis-eQTL at PCOS lead SNPs\n5  PCOS polygenic risk score [pending]")
arr(28,42,36,42); arr(50,36,50,31); arr(50,19,50,14)
arr(64,8,72,8); arr(64,8,72,25); arr(64,8,72,42)
 
fig.savefig(f"{F}/Figure1_design.png",bbox_inches="tight"); fig.savefig(f"{F}/Figure1_design.pdf",bbox_inches="tight"); plt.close()
 
# ---------- Fig 2 donors ----------
qc=pd.read_csv("/home/claude/geno/qc_samples.tsv",sep="\t",dtype={"IID":str})
fig=plt.figure(figsize=(7.6,4.4)); gsp=fig.add_gridspec(2,4,wspace=.8,hspace=.65,width_ratios=[1,1,1,1.35])
for i,(c,l,lgs) in enumerate([("amh","AMH (ng/mL)",False),("testosterone","Testosterone (pmol/L)",False),("fai","Free androgen index",True),("fg_score","Ferriman–Gallwey score",False),("bmi","BMI (kg/m²)",False),("homa_ir","HOMA-IR",True)]):
    ax=fig.add_subplot(gsp[i//3, i%3] if i<3 else gsp[1,i-3]); strip(ax,c,l,lgs); lab(ax,"ABCDEF"[i])
ax=fig.add_subplot(gsp[:,3])
don=dict(zip(key.array_iid,key.group))
qc=qc[~qc.IID.isin(__import__("os").environ.get("RELATED_PAIR_SAMPLE_IDS","").split(","))]  # a first-degree related pair of samples is omitted from the PCA panel; IDs held by the authors (comma-separated, set in the environment)
oth=qc[~qc.IID.isin(don)]
ax.scatter(oth.PC2,oth.PC3,s=10,color="#c4c7cc",lw=0,label="Genotyped cohort")
for g in ["control","PMOS"]:
    d=qc[qc.IID.map(don)==g]; ax.scatter(d.PC2,d.PC3,s=28,marker=M[g],color=C[g],edgecolor="white",lw=.6,label=f"iPSC donor, {g}")
ax.set_xlabel("Genotype PC2"); ax.set_ylabel("Genotype PC3"); ax.legend(frameon=False,fontsize=6.3,loc="upper center",bbox_to_anchor=(0.5,-0.13)); lab(ax,"G")
fig.savefig(f"{F}/Figure2_donors.png",bbox_inches="tight"); fig.savefig(f"{F}/Figure2_donors.pdf",bbox_inches="tight"); plt.close()
 
# ---------- Fig 3 transcriptome ----------
pcs=pd.read_csv(f"{O}/tables/pcs.tsv",sep="\t",index_col=0).loc[lg.columns]
rho=pd.read_csv(f"{O}/tables/pc_trait_rho.tsv",sep="\t",index_col=0); pp=pd.read_csv(f"{O}/tables/pc_trait_p.tsv",sep="\t",index_col=0)
perm=pd.read_csv(f"{O}/tables/S5b_global_permutation.tsv",sep="\t")
de=pd.read_csv(f"{O}/tables/S4_DE_PMOS_vs_control.tsv",sep="\t")
fig=plt.figure(figsize=(7.2,5.9)); gsp=fig.add_gridspec(2,2,wspace=.45,hspace=.62)
ax=fig.add_subplot(gsp[0,0])
for g in ["control","PMOS"]:
    d=pcs[key.group==g]; ax.scatter(d.PC1,d.PC2,s=34,marker=M[g],color=C[g],edgecolor="white",lw=.6,label="Control" if g=="control" else "PMOS")
for l,r in pcs.iterrows(): ax.annotate(l.replace("HEL",""),(r.PC1,r.PC2),fontsize=5.5,color="#555",xytext=(3,2),textcoords="offset points")
vlab=rho.index
ax.set_xlabel(f"PC1 ({vlab[0].split('(')[1][:-1]})"); ax.set_ylabel(f"PC2 ({vlab[1].split('(')[1][:-1]})"); ax.legend(frameon=False,fontsize=7); lab(ax,"A")
ax=fig.add_subplot(gsp[0,1])
im=ax.imshow(rho.values,cmap="RdBu_r",vmin=-.8,vmax=.8,aspect="auto")
ax.set_xticks(range(rho.shape[1])); ax.set_xticklabels(rho.columns,rotation=45,ha="right",fontsize=6.5)
ax.set_yticks(range(rho.shape[0])); ax.set_yticklabels(rho.index,fontsize=6.5)
for i in range(rho.shape[0]):
    for j in range(rho.shape[1]):
        if pp.values[i,j]<0.05: ax.text(j,i,"*",ha="center",va="center",fontsize=9,color=INK)
cb=fig.colorbar(im,ax=ax,fraction=.04,pad=.02); cb.set_label("Spearman ρ",fontsize=7); cb.ax.tick_params(labelsize=6)
ax.spines[:].set_visible(False); lab(ax,"B")
ax=fig.add_subplot(gsp[1,0])
P=perm.copy(); P=P.iloc[::-1]
y=np.arange(len(P))
ax.hlines(y,0,P.null_95pct,color="#dde0e4",lw=6,zorder=1,label="Permutation null, 0–95th pct")
ax.scatter(P.null_median,y,marker="|",s=60,color=MUTED,zorder=2,label="Null median")
ax.scatter(P.genes_p05,y,s=30,color=INK,zorder=3,label="Observed")
for yi,(o,p) in enumerate(zip(P.genes_p05,P.perm_p)): ax.text(max(o,P.null_95pct.iloc[yi])+60,yi,f"p = {p:.2f}",va="center",fontsize=6.5,color="#444")
ax.set_yticks(y); ax.set_yticklabels(P.trait); ax.set_xlabel("Genes associated at nominal p < 0.05"); ax.set_xlim(0,2800)
ax.legend(frameon=False,fontsize=6,loc="lower center",bbox_to_anchor=(0.45,-0.42),ncol=3); lab(ax,"C")
ax=fig.add_subplot(gsp[1,1])
ax.hist(de.p,bins=40,color=C["PMOS"],alpha=.85,edgecolor="white",lw=.4,label="PMOS vs control")
ax.hist(de.p_adjXIST,bins=40,histtype="step",color=INK,lw=1,label="+ XIST covariate")
ax.axhline(len(de)/40,color=MUTED,ls="--",lw=.8,label="Uniform (null) expectation")
ax.set_xlabel("Nominal p value"); ax.set_ylabel("Genes"); ax.set_ylim(0,620); ax.legend(frameon=False,fontsize=6.3,loc="upper center",ncol=2,bbox_to_anchor=(0.5,1.02)); lab(ax,"D")
fig.savefig(f"{F}/Figure3_transcriptome.png",bbox_inches="tight"); fig.savefig(f"{F}/Figure3_transcriptome.pdf",bbox_inches="tight"); plt.close()
 
# ---------- Fig 4 XIST ----------
x=np.log2(cnt.loc["XIST",lg.columns]+1)
order=x.sort_values().index
fig,axs=plt.subplots(1,3,figsize=(7.2,2.6),gridspec_kw=dict(width_ratios=[2.2,1,1.3],wspace=.55))
ax=axs[0]
ax.bar(range(len(order)),x[order],color=[C[key.loc[l,"group"]] for l in order],width=.7)
ax.set_xticks(range(len(order))); ax.set_xticklabels([l.replace("HEL","") for l in order],rotation=60,ha="right",fontsize=6.3)
ax.set_ylabel("XIST, log2(normalised count + 1)"); ax.axhline(8,color=MUTED,ls=":",lw=.8)
from matplotlib.patches import Patch
ax.legend(handles=[Patch(color=C["control"],label="Control"),Patch(color=C["PMOS"],label="PMOS")],frameon=False,fontsize=6.5,loc="upper left"); lab(ax,"A")
ax=axs[1]
for g in ["control","PMOS"]:
    v=x[key.group==g]; i=0 if g=="control" else 1
    ax.scatter(i+np.random.default_rng(1).uniform(-.12,.12,len(v)),v,s=24,marker=M[g],color=C[g],edgecolor="white",lw=.6)
    ax.hlines(v.median(),i-.25,i+.25,color=INK,lw=1.2)
p=stats.mannwhitneyu(x[key.group=="PMOS"],x[key.group=="control"]).pvalue
ax.set_xticks([0,1]); ax.set_xticklabels(["Control","PMOS"]); ax.set_xlim(-.6,1.6); ax.set_ylabel("XIST (log2)"); ax.set_title(f"p = {p:.3f}",fontsize=7.5,fontweight="normal",color="#444"); lab(ax,"B")
ax=axs[2]
k=np.log2(cnt.loc["KDM5C",lg.columns]+1)
for g in ["control","PMOS"]:
    m=key.group==g; ax.scatter(x[m],k[m],s=24,marker=M[g],color=C[g],edgecolor="white",lw=.6)
r,pv=stats.spearmanr(x,k); ax.set_xlabel("XIST (log2)"); ax.set_ylabel("KDM5C (log2)"); ax.set_title(f"ρ = {r:.2f}, p = {pv:.3f}",fontsize=7.5,fontweight="normal",color="#444"); lab(ax,"C")
fig.savefig(f"{F}/Figure4_XIST.png",bbox_inches="tight"); fig.savefig(f"{F}/Figure4_XIST.pdf",bbox_inches="tight"); plt.close()
 
# ---------- Fig 5 GWAS genes ----------
G=pd.read_csv(f"{O}/tables/S7_PCOS_GWAS_genes_expression.tsv",sep="\t")
cat={"Endocrine, gonadal & germline":["FSHR","FSHB","AMH","SHBG","INHBB","GATA4","CYP3A7","CYP3A5","HHAT","TOPAZ1","ZBTB16","TEX41"],
     "DNA repair / genome maintenance":["MSH6","CHEK2","FANCC","RAD50","NEIL2","HELB"],
     "Growth, signalling & other":["THADA","DENND1A","YAP1","HMGA2","ERBB3","ERBB4","MYO10","TOX3","FTO","MAF","KRR1","ID4","PRSS23","IRF1","MAPRE1"]}
cc={"Endocrine, gonadal & germline":"#e87ba4","DNA repair / genome maintenance":"#1baf7a","Growth, signalling & other":"#2a78d6"}
G["category"]=G.gene.map({g:k for k,v in cat.items() for g in v})
G["log2mean"]=np.log2(G.mean_norm_count.fillna(0)+1)
G=G.sort_values(["category","log2mean"],ascending=[True,True])
fig=plt.figure(figsize=(7.2,5.2)); gsp=fig.add_gridspec(1,2,width_ratios=[1.6,1],wspace=.45)
ax=fig.add_subplot(gsp[0,0])
y=np.arange(len(G)); allm=np.log2(cnt.loc[lg.index,lg.columns].mean(1)+1)
q=np.percentile(allm,[25,50,75])
ax.axvspan(q[0],q[2],color="#eef0f3",zorder=0); ax.axvline(q[1],color=MUTED,lw=.8,ls="--",zorder=1)
ax.axvline(np.log2(11),color="#c4c7cc",lw=.8,zorder=1)
ax.hlines(y,0,G.log2mean,color="#d0d3d8",lw=1,zorder=2)
ax.scatter(G.log2mean,y,s=26,color=[cc[c] for c in G.category],zorder=3)
ax.set_yticks(y); ax.set_yticklabels(G.gene,fontsize=6.3,style="italic"); ax.set_xlabel("Mean expression, log2(normalised count + 1)")
ax.text(q[1],len(G)-.2,"expressed-gene\nmedian / IQR",fontsize=5.8,color=MUTED,ha="center",va="bottom")
ax.text(np.log2(11),len(G)-.2,"expression\nthreshold",fontsize=5.8,color=MUTED,ha="center",va="bottom")
ax.legend(handles=[Patch(color=v,label=k) for k,v in cc.items()],frameon=False,fontsize=6.3,loc="upper center",bbox_to_anchor=(0.45,-0.09),ncol=2); lab(ax,"A")
ax=fig.add_subplot(gsp[0,1])
sc=pd.read_csv(f"{O}/tables/geneset_score.tsv",sep="\t",index_col=0).loc[lg.columns]
gt=pd.read_csv(f"{O}/tables/S7b_geneset_tests.tsv",sep="\t")
for i,g in enumerate(["control","PMOS"]):
    v=sc.score[sc.group==g]; ax.scatter(i+np.random.default_rng(2).uniform(-.12,.12,len(v)),v,s=26,marker=M[g],color=C[g],edgecolor="white",lw=.6); ax.hlines(v.median(),i-.25,i+.25,color=INK,lw=1.2)
ax.set_xticks([0,1]); ax.set_xticklabels(["Control","PMOS"]); ax.set_xlim(-.6,1.6); ax.set_ylabel("PCOS GWAS gene-set score (mean z)")
pr=gt.set_index("test").loc["PCOS GWAS gene-set score vs PMOS status","p_perm_matched"]
ax.set_title(f"24 expressed genes\nexpression-matched permutation p = {pr:.2f}",fontsize=7,fontweight="normal",color="#444"); lab(ax,"B")
fig.savefig(f"{F}/Figure5_GWAS_genes.png",bbox_inches="tight"); fig.savefig(f"{F}/Figure5_GWAS_genes.pdf",bbox_inches="tight"); plt.close()
 
# ---------- Fig 6 genetics ----------
gd=pd.read_csv("/home/claude/geno/donor_genotypes_6snps.tsv",sep="\t",index_col=0).loc[lg.columns]
eq=pd.read_csv("/home/claude/geno/targeted_eqtl_6snps.tsv",sep="\t")
fig=plt.figure(figsize=(7.2,5.0)); gsp=fig.add_gridspec(2,3,wspace=.55,hspace=.6)
ax=fig.add_subplot(gsp[0,0])
vals=[29,7,5]; labs=["Lead variants\n(Moolhuijsen 2026)","On GSA array","Passed QC"]
ax.barh([2,1,0],vals,color=["#c4c7cc","#9aa0a6",INK],height=.55)
for yy,v in zip([2,1,0],vals): ax.text(v+.6,yy,str(v),va="center",fontsize=7.5)
ax.set_yticks([2,1,0]); ax.set_yticklabels(labs,fontsize=6.5); ax.set_xlim(0,34); ax.set_xlabel("Variants"); lab(ax,"A")
def eqp(ax,snp,gene,letter):
    d=gd[snp]; e=np.log2(cnt.loc[gene,lg.columns]+1); ok=d.notna()
    rng=np.random.default_rng(5)
    for g in ["control","PMOS"]:
        m=ok&(key.group==g); ax.scatter(d[m]+rng.uniform(-.1,.1,m.sum()),e[m],s=24,marker=M[g],color=C[g],edgecolor="white",lw=.6)
    for k in (0,1,2):
        v=e[ok&(d==k)]
        if len(v): ax.hlines(v.median(),k-.25,k+.25,color=INK,lw=1.2)
    r=eq[(eq.snp==snp)&(eq.gene==gene)].iloc[0]
    ax.set_xticks([0,1,2]); ax.set_xlabel(f"{snp} risk alleles"); ax.set_ylabel(f"{gene} (log2)",style="italic")
    ax.set_title(f"β = {r.beta_log2_per_risk:+.2f}, p = {r.p:.3f}\n(q = {r.q_BH:.2f})",fontsize=7,fontweight="normal",color="#444"); lab(ax,letter)
eqp(fig.add_subplot(gsp[0,1]),"rs10876864","ERBB3","B"); eqp(fig.add_subplot(gsp[0,2]),"rs10876864","RPS26","C")
ax=fig.add_subplot(gsp[1,:2])
E=eq[eq.status=="tested"].copy(); E["lab"]=E.gene+" · "+E.snp
E=E.sort_values("beta_log2_per_risk")
yy=np.arange(len(E))
ax.errorbar(E.beta_log2_per_risk,yy,xerr=1.96*E.se,fmt="o",color=INK,ms=4,elinewidth=.8,capsize=0)
ax.axvline(0,color=MUTED,lw=.8)
ax.set_yticks(yy); ax.set_yticklabels(E.lab,fontsize=6.3); ax.set_xlabel("log2 expression change per PCOS risk allele (95% CI)"); lab(ax,"D")
ax=fig.add_subplot(gsp[1,2]); ax.set_xticks([]); ax.set_yticks([])
for s in ax.spines.values(): s.set_visible(True); s.set_linestyle("--"); s.set_color("#9aa0a6")
ax.text(.5,.5,"PCOS polygenic risk score\n(PRS-CS, ~1.1 M variants)\nacross 140 genotyped women\n\n[PENDING: imputation and\nGWAS weights]",ha="center",va="center",fontsize=7,color="#666",transform=ax.transAxes); lab(ax,"E")
fig.savefig(f"{F}/Figure6_genetics.png",bbox_inches="tight"); fig.savefig(f"{F}/Figure6_genetics.pdf",bbox_inches="tight"); plt.close()
print("done")
