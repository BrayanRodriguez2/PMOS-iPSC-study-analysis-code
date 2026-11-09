"""v3 figures (2 Nov 2026): Figures 3-5 and Supplementary Figures S1-S3.
Usage: python3 figs_v3.py <package_dir> <counts.tabular> <an_restore_out_dir>"""
import sys
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from scipy import stats
from scipy.cluster.hierarchy import linkage, leaves_list
PKG, COUNTS, A = sys.argv[1:4]
sys.path.insert(0, PKG + "/5_Pipeline_CSC/scripts")
from common import read_counts
R = PKG + "/4_Analysis_results"; F = PKG + "/2_Figures"
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":8,"axes.spines.top":False,"axes.spines.right":False,
  "axes.linewidth":0.6,"xtick.major.width":0.6,"ytick.major.width":0.6,"axes.edgecolor":"#555","axes.labelcolor":"#222",
  "xtick.color":"#444","ytick.color":"#444","axes.titlesize":9,"axes.titleweight":"bold","savefig.dpi":300,"pdf.fonttype":42})
C={"PMOS":"#eb6834","control":"#2a78d6"}; INK="#222"; MUTED="#888"
BM={0:"o",1:"^"}  # lean circle, obese triangle
def lab(ax,s,x=-0.16): ax.text(x,1.06,s,transform=ax.transAxes,fontsize=11,fontweight="bold",va="bottom")
def save(fig,name): fig.savefig(f"{F}/{name}.png",bbox_inches="tight"); fig.savefig(f"{F}/{name}.pdf",bbox_inches="tight"); plt.close(fig)
 
key=pd.read_csv(f"{R}/metadata/sample_key_14_ipsc_donors.csv").set_index("hel_line")
cnt=read_counts(COUNTS); la=pd.read_csv(f"{A}/line_annot.tsv",sep="\t").set_index("line")
lines=list(la.index); cnt=cnt[lines]
lg=np.log2(cnt[(cnt>=10).sum(axis=1)>=7]+1)
pcs=pd.read_csv(f"{R}/transcriptome/pcs.tsv",sep="\t",index_col=0).loc[lines]
perm=pd.read_csv(f"{R}/transcriptome/S5b_global_permutation.tsv",sep="\t")
vol=pd.read_csv(f"{A}/volcano.tsv.gz",sep="\t")
Cm=pd.read_csv(f"{A}/sample_correlation.tsv",sep="\t",index_col=0)
ct=pd.read_csv(f"{A}/sample_correlation_tests.tsv",sep="\t").set_index("grouping")
 
# ================= Figure 3 =================
fig=plt.figure(figsize=(7.4,7.0)); gsp=fig.add_gridspec(2,2,wspace=.42,hspace=.42,width_ratios=[1,1.05])
ax=fig.add_subplot(gsp[0,0])
for l,r in pcs.iterrows():
    g=la.loc[l,"group"]; ax.scatter(r.PC1,r.PC2,s=36,marker=BM[int(la.loc[l,"obese"])],color=C[g],edgecolor="white",lw=.6,zorder=3)
OFF={"HEL351.3":(-19,3),"HEL364.2":(-4,5),"HEL350.3":(4,-5),"HEL359.1":(-6,5),"HEL361.2":(4,-7),"HEL363.2":(4,-4)}
for l,r in pcs.iterrows(): ax.annotate(l.replace("HEL",""),(r.PC1,r.PC2),fontsize=5.5,color="#555",xytext=OFF.get(l,(3,2)),textcoords="offset points")
ax.set_xlim(-32,66); ax.set_xlabel("PC1 (43%)"); ax.set_ylabel("PC2 (9%)")
h=[Line2D([],[],ls="",marker="s",color=C["control"],label="Control"),Line2D([],[],ls="",marker="s",color=C["PMOS"],label="PMOS"),
   Line2D([],[],ls="",marker="o",color="#777",label="BMI <25"),Line2D([],[],ls="",marker="^",color="#777",label="BMI >30")]
ax.legend(handles=h,frameon=False,fontsize=6.3,loc="upper right",ncol=2,handletextpad=.2,columnspacing=.8); lab(ax,"A")
# B correlation heatmap, clustered
ax=fig.add_subplot(gsp[0,1])
order=leaves_list(linkage(1-Cm.values[np.triu_indices(14,1)],"average")); o=[lines[i] for i in order]
M=Cm.loc[o,o].values.copy(); np.fill_diagonal(M,np.nan)
im=ax.imshow(M,cmap="viridis",vmin=np.nanmin(M),vmax=np.nanmax(M),aspect="equal",extent=(0,14,14,0))
ax.set_xticks(np.arange(14)+.5); ax.set_xticklabels([l.replace("HEL","") for l in o],rotation=90,fontsize=5.3)
ax.set_yticks(np.arange(14)+.5); ax.set_yticklabels([l.replace("HEL","") for l in o],fontsize=5.3)
for j,l in enumerate(o):
    ax.add_patch(plt.Rectangle((j,-1.0),1,.8,color=C[la.loc[l,"group"]],clip_on=False))
    ax.add_patch(plt.Rectangle((j,-1.9),1,.8,color="#9aa0a6" if la.loc[l,"obese"] else "#e3e5e8",clip_on=False))
    ax.add_patch(plt.Rectangle((j,-2.8),1,.8,color="#6b4fa0" if la.loc[l,"xist_log2"]<8 else "#e3e5e8",clip_on=False))
for yy,t in [(-0.6,"Group"),(-1.5,"BMI >30"),(-2.4,"XIST eroded")]: ax.text(14.3,yy,t,fontsize=5.6,va="center",color="#444")
ax.set_ylim(14,-2.8); ax.spines[:].set_visible(False); ax.tick_params(length=0)
cb=fig.colorbar(im,ax=ax,fraction=.04,pad=.25); cb.set_label("Pearson r",fontsize=7); cb.ax.tick_params(labelsize=6)
lab(ax,"B",-.12)
# C volcano
ax=fig.add_subplot(gsp[1,0])
v=vol.copy(); v["nl"]=-np.log10(v.p); up=(v.p<.05)&(v.log2FC>0); dn=(v.p<.05)&(v.log2FC<0)
ax.scatter(v.log2FC[~(up|dn)],v.nl[~(up|dn)],s=3,color="#cfd2d6",lw=0,rasterized=True)
ax.scatter(v.log2FC[up],v.nl[up],s=5,color=C["PMOS"],lw=0,rasterized=True,label=f"higher in PMOS (p<0.05, n={up.sum()})")
ax.scatter(v.log2FC[dn],v.nl[dn],s=5,color=C["control"],lw=0,rasterized=True,label=f"lower in PMOS (p<0.05, n={dn.sum()})")
ax.axhline(-np.log10(.05),color=MUTED,ls="--",lw=.7)
VO={"C11orf54":(0,6,"center"),"NPHP3":(-3,-8,"right")}
for _,r in v.nsmallest(8,"p").iterrows():
    neg=r.log2FC<0; dx,dy,ha=VO.get(r.gene,(-3 if neg else 3,1,"right" if neg else "left"))
    ax.annotate(r.gene,(r.log2FC,r.nl),fontsize=5.5,style="italic",ha=ha,xytext=(dx,dy),textcoords="offset points")
ax.set_xlabel("log2 fold change, PMOS vs control"); ax.set_ylabel("−log10 p"); ax.legend(frameon=False,fontsize=5.8,loc="lower left",bbox_to_anchor=(0,-0.42),ncol=1)
ax.set_title("16,489 genes; none at FDR < 0.1",fontsize=7,fontweight="normal",color="#444"); lab(ax,"C")
# D permutation null
ax=fig.add_subplot(gsp[1,1])
P=perm.iloc[::-1]; y=np.arange(len(P))
ax.hlines(y,0,P.null_95pct,color="#dde0e4",lw=6,zorder=1,label="Null, 0–95th pct")
ax.scatter(P.null_median,y,marker="|",s=60,color=MUTED,zorder=2,label="Null median")
ax.scatter(P.genes_p05,y,s=28,color=INK,zorder=3,label="Observed")
for yi,(ob,p) in enumerate(zip(P.genes_p05,P.perm_p)): ax.text(max(ob,P.null_95pct.iloc[yi])+60,yi,f"p = {p:.2f}",va="center",fontsize=6.3,color="#444")
ax.set_yticks(y); ax.set_yticklabels(P.trait); ax.set_xlabel("Genes associated at nominal p < 0.05"); ax.set_xlim(0,2800)
ax.legend(frameon=False,fontsize=5.8,loc="lower center",bbox_to_anchor=(0.45,-0.4),ncol=3); lab(ax,"D",-.3)
save(fig,"Figure3_transcriptome")
 
# ================= Figure 4 =================
vpt=pd.read_csv(f"{A}/variance_partition_per_gene.tsv.gz",sep="\t",index_col=0)
vps=pd.read_csv(f"{A}/variance_partition_summary.tsv",sep="\t").set_index("factor")
bs=pd.read_csv(f"{A}/bmi_adjusted_summary.tsv",sep="\t")
fig=plt.figure(figsize=(7.4,6.0)); gsp=fig.add_gridspec(2,2,wspace=.38,hspace=.55,width_ratios=[1.45,1])
ax=fig.add_subplot(gsp[0,0])
cols=["PMOS status","BMI","Age","XIST (log2)","Residual"]; vc=["#eb6834","#e0a33a","#1baf7a","#6b4fa0","#9aa0a6"]
parts=ax.violinplot([100*vpt[c].values for c in cols],positions=range(5),widths=.8,showextrema=False)
for b,cl in zip(parts["bodies"],vc): b.set_facecolor(cl); b.set_alpha(.45); b.set_edgecolor("none")
for i,c in enumerate(cols):
    m=vps.loc[c,"median_pct"]; ax.hlines(m,i-.3,i+.3,color=INK,lw=1.4,zorder=3)
    ax.text(i,max(m,0)+(4 if c!="Residual" else -14),f"{m:.1f}%",ha="center",fontsize=6.3)
    if c!="Residual":
        ax.vlines(i+.36,0,vps.loc[c,"null_95pct"],color=MUTED,lw=2.2,alpha=.6)
        ax.scatter(i+.36,vps.loc[c,"null_median_pct"],marker="_",s=40,color=INK,zorder=4)
        ax.text(i,-13,f"p = {vps.loc[c,'perm_p']:.2f}",ha="center",fontsize=6,color="#444")
ax.set_xticks(range(5)); ax.set_xticklabels(["PMOS","BMI","Age","XIST","Residual"]); ax.set_ylim(-17,100)
ax.set_ylabel("Variance explained per gene (%)"); ax.axhline(0,color="#bbb",lw=.5)
ax.legend(handles=[Line2D([],[],color=INK,lw=1.4,label="Observed median"),Line2D([],[],color=MUTED,lw=2.2,alpha=.6,label="Permutation null, 0–95th pct (median tick)")],
          frameon=False,fontsize=5.8,loc="upper left"); lab(ax,"A",-.12)
ax=fig.add_subplot(gsp[0,1])
u=int(bs.genes_p05[0]); a=int(bs.genes_p05[1]); ov=int(bs.overlap_with_main[1]); ex=int(bs.expected_by_chance[0])
bars=[("Unadjusted",u,C["PMOS"]),("BMI-\nadjusted",a,"#e0a33a"),("Overlap",ov,"#9aa0a6"),("Expected\nby chance",ex,"#e3e5e8")]
for i,(t,vv,cl) in enumerate(bars): ax.bar(i,vv,color=cl,width=.7,edgecolor="#999" if i==3 else "none",lw=.5); ax.text(i,vv+15,str(vv),ha="center",fontsize=6.5)
ax.set_xticks(range(4)); ax.set_xticklabels([b[0] for b in bars],fontsize=6.5); ax.set_ylabel("Genes, PMOS term p < 0.05"); ax.set_ylim(0,950)
ax.set_title(f"{100*ov/u:.0f}% of unadjusted hits retained; none at FDR < 0.1",fontsize=6.6,fontweight="normal",color="#444"); lab(ax,"B")
x=np.log2(cnt.loc["XIST"]+1); order=x.sort_values().index
ax=fig.add_subplot(gsp[1,0])
ax.bar(range(14),x[order],color=[C[la.loc[l,"group"]] for l in order],width=.7)
ax.set_xticks(range(14)); ax.set_xticklabels([l.replace("HEL","") for l in order],rotation=60,ha="right",fontsize=6.2)
ax.set_ylabel("XIST, log2(normalised count + 1)"); ax.axhline(8,color=MUTED,ls=":",lw=.8)
ax.legend(handles=[Patch(color=C["control"],label="Control"),Patch(color=C["PMOS"],label="PMOS")],frameon=False,fontsize=6.3,loc="upper left")
ax.set_title("3/7 control vs 0/7 PMOS lines eroded (Fisher p = 0.19)",fontsize=6.6,fontweight="normal",color="#444"); lab(ax,"C",-.12)
ax=fig.add_subplot(gsp[1,1]); k=np.log2(cnt.loc["KDM5C"]+1)
for g in ["control","PMOS"]:
    m=la.group==g; ax.scatter(x[m.values],k[m.values],s=24,marker="s" if g=="control" else "o",color=C[g],edgecolor="white",lw=.6)
r,pv=stats.spearmanr(x,k); ax.set_xlabel("XIST (log2)"); ax.set_ylabel("KDM5C (log2)")
ax.set_title(f"ρ = {r:.2f}, p = {pv:.3f}",fontsize=7,fontweight="normal",color="#444"); lab(ax,"D")
save(fig,"Figure4_variation_XIST")
 
# ================= Figure 5 =================
G=pd.read_csv(f"{R}/transcriptome/S7_PCOS_GWAS_genes_expression.tsv",sep="\t")
pl=pd.read_csv(f"{A}/protein_pcos_locus_genes.tsv",sep="\t").set_index("gene")
cat={"Endocrine, gonadal & germline":["FSHR","FSHB","AMH","SHBG","INHBB","GATA4","CYP3A7","CYP3A5","HHAT","TOPAZ1","ZBTB16","TEX41"],
     "DNA repair / genome maintenance":["MSH6","CHEK2","FANCC","RAD50","NEIL2","HELB"],
     "Growth, signalling & other":["THADA","DENND1A","YAP1","HMGA2","ERBB3","ERBB4","MYO10","TOX3","FTO","MAF","KRR1","ID4","PRSS23","IRF1","MAPRE1"]}
cc={"Endocrine, gonadal & germline":"#e87ba4","DNA repair / genome maintenance":"#1baf7a","Growth, signalling & other":"#2a78d6"}
G["category"]=G.gene.map({g:k for k,v in cat.items() for g in v}); G["log2mean"]=np.log2(G.mean_norm_count.fillna(0)+1)
G=G.sort_values(["category","log2mean"])
fig=plt.figure(figsize=(7.2,5.2)); gsp=fig.add_gridspec(1,2,width_ratios=[1.6,1],wspace=.45)
ax=fig.add_subplot(gsp[0,0]); y=np.arange(len(G)); allm=np.log2(cnt.loc[lg.index].mean(axis=1)+1); q=np.percentile(allm,[25,50,75])
ax.axvspan(q[0],q[2],color="#eef0f3",zorder=0); ax.axvline(q[1],color=MUTED,lw=.8,ls="--",zorder=1); ax.axvline(np.log2(11),color="#c4c7cc",lw=.8,zorder=1)
ax.hlines(y,0,G.log2mean,color="#d0d3d8",lw=1,zorder=2); ax.scatter(G.log2mean,y,s=26,color=[cc[c] for c in G.category],zorder=3)
px=14.2
for yi,gname in zip(y,G.gene):
    if bool(pl.loc[gname,"protein_detected_HipSci"]): ax.scatter(px,yi,marker="s",s=16,color=INK,clip_on=False)
ax.text(px,len(G)-.2,"protein\n(HipSci)",fontsize=5.8,color=INK,ha="center",va="bottom")
ax.set_xlim(0,14.8); ax.set_yticks(y); ax.set_yticklabels(G.gene,fontsize=6.3,style="italic"); ax.set_xlabel("Mean expression, log2(normalised count + 1)")
ax.text(q[1],len(G)-.2,"expressed-gene\nmedian / IQR",fontsize=5.8,color=MUTED,ha="center",va="bottom")
ax.text(np.log2(11),len(G)-.2,"expression\nthreshold",fontsize=5.8,color=MUTED,ha="center",va="bottom")
ax.legend(handles=[Patch(color=v,label=k) for k,v in cc.items()]+[Line2D([],[],ls="",marker="s",color=INK,label="protein detected in HipSci iPSCs")],
          frameon=False,fontsize=6.2,loc="upper center",bbox_to_anchor=(0.45,-0.09),ncol=2); lab(ax,"A")
ax=fig.add_subplot(gsp[0,1])
sc=pd.read_csv(f"{R}/transcriptome/geneset_score.tsv",sep="\t",index_col=0); gt=pd.read_csv(f"{R}/transcriptome/S7b_geneset_tests.tsv",sep="\t")
for i,g in enumerate(["control","PMOS"]):
    vv=sc.score[sc.group==g]; ax.scatter(i+np.random.default_rng(2).uniform(-.12,.12,len(vv)),vv,s=26,marker="s" if g=="control" else "o",color=C[g],edgecolor="white",lw=.6); ax.hlines(vv.median(),i-.25,i+.25,color=INK,lw=1.2)
ax.set_xticks([0,1]); ax.set_xticklabels(["Control","PMOS"]); ax.set_xlim(-.6,1.6); ax.set_ylabel("PCOS GWAS gene-set score (mean z)")
pr=gt.set_index("test").loc["PCOS GWAS gene-set score vs PMOS status","p_perm_matched"]
ax.set_title(f"24 expressed genes\nexpression-matched permutation p = {pr:.2f}",fontsize=7,fontweight="normal",color="#444"); lab(ax,"B")
save(fig,"Figure5_GWAS_genes")
 
# ================= Supplementary Figure S1 =================
rho=pd.read_csv(f"{R}/transcriptome/pc_trait_rho.tsv",sep="\t",index_col=0); pp=pd.read_csv(f"{R}/transcriptome/pc_trait_p.tsv",sep="\t",index_col=0)
de=pd.read_csv(f"{R}/transcriptome/S4_DE_PMOS_vs_control.tsv",sep="\t"); pa=pd.read_csv(f"{A}/pcs_allgenes.tsv",sep="\t",index_col=0).loc[lines]
fig=plt.figure(figsize=(7.8,2.9)); gsp=fig.add_gridspec(1,3,wspace=.85,width_ratios=[1,1.15,1])
ax=fig.add_subplot(gsp[0,0])
for l,r in pa.iterrows(): ax.scatter(r.PC1,r.PC2,s=30,marker=BM[int(la.loc[l,"obese"])],color=C[la.loc[l,"group"]],edgecolor="white",lw=.6)
OA={"HEL364.2":(4,3),"HEL350.3":(4,-7),"HEL358.1":(-3,-9),"HEL363.2":(4,2),"HEL351.3":(4,2)}
for l,r in pa.iterrows(): ax.annotate(l.replace("HEL",""),(r.PC1,r.PC2),fontsize=5,color="#555",xytext=OA.get(l,(3,2)),textcoords="offset points")
ax.set_xlabel("PC1 (35%)"); ax.set_ylabel("PC2 (10%)"); ax.set_title("All 16,489 expressed genes",fontsize=7,fontweight="normal",color="#444"); lab(ax,"A",-.25)
ax=fig.add_subplot(gsp[0,1]); im=ax.imshow(rho.values,cmap="RdBu_r",vmin=-.8,vmax=.8,aspect="auto")
ax.set_xticks(range(rho.shape[1])); ax.set_xticklabels(rho.columns,rotation=45,ha="right",fontsize=6)
ax.set_yticks(range(rho.shape[0])); ax.set_yticklabels(rho.index,fontsize=6)
for i in range(rho.shape[0]):
    for j in range(rho.shape[1]):
        if pp.values[i,j]<0.05: ax.text(j,i,"*",ha="center",va="center",fontsize=9,color=INK)
cb=fig.colorbar(im,ax=ax,fraction=.05,pad=.03); cb.set_label("Spearman ρ",fontsize=6.5); cb.ax.tick_params(labelsize=5.5); ax.spines[:].set_visible(False); lab(ax,"B",-.35)
ax=fig.add_subplot(gsp[0,2])
ax.hist(de.p,bins=40,color=C["PMOS"],alpha=.85,edgecolor="white",lw=.4,label="PMOS vs control")
ax.hist(de.p_adjXIST,bins=40,histtype="step",color=INK,lw=1,label="+ XIST covariate")
ax.axhline(len(de)/40,color=MUTED,ls="--",lw=.8,label="Uniform expectation")
ax.set_xlabel("Nominal p value"); ax.set_ylabel("Genes"); ax.set_ylim(0,720); ax.legend(frameon=False,fontsize=5.6,loc="upper center",ncol=1,bbox_to_anchor=(0.5,1.12)); lab(ax,"C",-.3)
save(fig,"FigureS1_PCA_traits_pvalues")
 
# ================= Supplementary Figure S2: top-40 heatmap =================
t40=pd.read_csv(f"{A}/top40_genes.tsv",sep="\t"); X=lg.loc[t40.gene]
ordl=[l for l in lines if la.loc[l,"group"]=="control"]+[l for l in lines if la.loc[l,"group"]=="PMOS"]
Zs=X[ordl].sub(X[ordl].mean(axis=1),axis=0).div(X[ordl].std(axis=1),axis=0)
fig,ax=plt.subplots(figsize=(5.2,7.2)); im=ax.imshow(Zs.values,cmap="RdBu_r",vmin=-2.2,vmax=2.2,aspect="auto")
ax.set_yticks(range(40)); ax.set_yticklabels([f"{g}  (p={p:.0e})" for g,p in zip(t40.gene,t40.p)],fontsize=5.4)
ax.set_xticks(range(14)); ax.set_xticklabels([l.replace("HEL","") for l in ordl],rotation=90,fontsize=6)
for j,l in enumerate(ordl): ax.add_patch(plt.Rectangle((j-.5,-1.6),1,.9,color=C[la.loc[l,"group"]],clip_on=False))
ax.set_ylim(39.5,-1.7); ax.spines[:].set_visible(False); ax.tick_params(length=0)
cb=fig.colorbar(im,ax=ax,fraction=.04,pad=.02); cb.set_label("z-score, log2 expression",fontsize=7)
ax.set_title("Top 40 genes by nominal p (PMOS vs control); none at FDR < 0.1",fontsize=7,fontweight="normal",color="#444",pad=14)
save(fig,"FigureS2_top40_heatmap")
 
# ================= Supplementary Figure S3: protein =================
dec=pd.read_csv(f"{A}/protein_by_decile.tsv",sep="\t"); cand=pd.read_csv(f"{A}/protein_candidates.tsv",sep="\t")
fig=plt.figure(figsize=(7.4,5.6)); gsp=fig.add_gridspec(1,2,wspace=.55,width_ratios=[1,1.5])
ax=fig.add_subplot(gsp[0,0]); ax.bar(dec.expression_decile,dec.pct_detected,color="#2a78d6",width=.7)
for d_,p_ in zip(dec.expression_decile,dec.pct_detected): ax.text(d_,p_+1.5,f"{p_:.0f}",ha="center",fontsize=5.8)
ax.set_xlabel("Expression decile in our iPSCs (1 = lowest)"); ax.set_ylabel("Genes with protein quantified in HipSci (%)")
ax.set_xticks(range(1,11)); ax.set_title("5,510 of 16,489 expressed genes (33%)",fontsize=7,fontweight="normal",color="#444"); lab(ax,"A",-.25)
ax=fig.add_subplot(gsp[0,1]); X=lg.loc[cand.gene,ordl]; Zs=X.sub(X.mean(axis=1),axis=0).div(X.std(axis=1),axis=0)
im=ax.imshow(Zs.values,cmap="RdBu_r",vmin=-2.2,vmax=2.2,aspect="auto")
ax.set_yticks(range(len(cand))); ax.set_yticklabels([f"{g} (r={r:.2f})" for g,r in zip(cand.gene,cand.hipsci_rna_protein_corr)],fontsize=5.2,style="italic")
ax.set_xticks(range(14)); ax.set_xticklabels([l.replace("HEL","") for l in ordl],rotation=90,fontsize=5.8)
for j,l in enumerate(ordl): ax.add_patch(plt.Rectangle((j-.5,-1.5),1,.8,color=C[la.loc[l,"group"]],clip_on=False))
ax.set_ylim(len(cand)-.5,-1.6); ax.spines[:].set_visible(False); ax.tick_params(length=0)
cb=fig.colorbar(im,ax=ax,fraction=.04,pad=.02); cb.set_label("z-score",fontsize=6.5)
ax.set_title(f"{len(cand)} genes: nominal p < 0.05 and HipSci RNA–protein r ≥ 0.3",fontsize=7,fontweight="normal",color="#444",pad=12); lab(ax,"B",-.35)
save(fig,"FigureS3_protein_crossref")
print("done")
