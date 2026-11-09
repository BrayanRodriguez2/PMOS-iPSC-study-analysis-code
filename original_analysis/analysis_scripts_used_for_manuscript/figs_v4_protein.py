"""v4 (4 Nov 2026): regenerates the two figures that depend on the HipSci protein cross-reference,
Figure 5 (PCOS GWAS-locus genes; black squares = protein detected) and Supplementary Figure S3,
from the v4 results of an_protein_v4.py. All other figures are unchanged from v3 (figs_v3.py).
Figure 5 code is identical to figs_v3.py except that it reads the v4 locus table.
 
Usage: python3 figs_v4_protein.py <package_dir> <counts.tabular> <v4_results_dir> <figure_out_dir>
"""
import sys, os
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
PKG, COUNTS, V4, F = sys.argv[1:5]
os.makedirs(F, exist_ok=True)
sys.path.insert(0, PKG + "/5_Pipeline_CSC/scripts")
from common import read_counts
R = PKG + "/4_Analysis_results"
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":8,"axes.spines.top":False,"axes.spines.right":False,
  "axes.linewidth":0.6,"xtick.major.width":0.6,"ytick.major.width":0.6,"axes.edgecolor":"#555","axes.labelcolor":"#222",
  "xtick.color":"#444","ytick.color":"#444","axes.titlesize":9,"axes.titleweight":"bold","savefig.dpi":300,"pdf.fonttype":42})
C={"PMOS":"#eb6834","control":"#2a78d6"}; INK="#222"; MUTED="#888"
def lab(ax,s,x=-0.16): ax.text(x,1.06,s,transform=ax.transAxes,fontsize=11,fontweight="bold",va="bottom")
def save(fig,name): fig.savefig(f"{F}/{name}.png",bbox_inches="tight"); fig.savefig(f"{F}/{name}.pdf",bbox_inches="tight"); plt.close(fig)
 
key=pd.read_csv(f"{R}/metadata/sample_key_14_ipsc_donors.csv").set_index("hel_line")
lines=list(key.sort_values(["group","bmi"]).index)
cnt=read_counts(COUNTS)[lines]
lg=np.log2(cnt[(cnt>=10).sum(axis=1)>=7]+1)
la=pd.DataFrame({"group":key.loc[lines,"group"],"obese":(key.loc[lines,"bmi"]>30).astype(int)},index=lines)
ordl=[l for l in lines if la.loc[l,"group"]=="control"]+[l for l in lines if la.loc[l,"group"]=="PMOS"]
 
# ================= Figure 5 (unchanged code; v4 protein markers) =================
G=pd.read_csv(f"{R}/transcriptome/S7_PCOS_GWAS_genes_expression.tsv",sep="\t")
pl=pd.read_csv(f"{V4}/S13b_protein_PCOS_locus_genes.tsv",sep="\t").set_index("gene")
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
 
# ================= Supplementary Figure S3 (v4) =================
dec=pd.read_csv(f"{V4}/S13a_protein_by_expression_decile.tsv",sep="\t")
cand=pd.read_csv(f"{V4}/S13c_protein_supported_candidates.tsv",sep="\t")
go=pd.read_csv(f"{V4}/S13e_GO_enrichment.tsv",sep="\t").dropna(subset=["p_adj_gSCS"])
d=pd.read_csv(f"{V4}/S13d_protein_crossref_all_genes.tsv.gz",sep="\t")
pcg=d[d.protein_coding]; ndet=int(pcg.protein_detected_HipSci.sum()); npc=len(pcg)
 
fig=plt.figure(figsize=(7.6,8.6)); gsp=fig.add_gridspec(2,2,wspace=.6,hspace=.42,width_ratios=[1,1.25],height_ratios=[1,1.05])
# A: detection by decile
ax=fig.add_subplot(gsp[0,0])
ax.bar(dec.expression_decile,dec.pct_detected_of_protein_coding,color="#2a78d6",width=.7,label="of protein-coding genes")
ax.scatter(dec.expression_decile,dec.pct_detected,marker="_",s=70,color=INK,zorder=3,label="of all expressed genes")
for d_,p_ in zip(dec.expression_decile,dec.pct_detected_of_protein_coding): ax.text(d_,p_+2,f"{p_:.0f}",ha="center",fontsize=5.8)
ax.set_xlabel("Expression decile in our iPSCs (1 = lowest)"); ax.set_ylabel("Genes with protein detected in HipSci (%)")
ax.set_xticks(range(1,11)); ax.set_ylim(0,108)
ax.legend(frameon=False,fontsize=5.8,loc="upper left")
ax.set_title(f"{ndet:,} of {npc:,} expressed protein-coding genes ({100*ndet/npc:.0f}%)",fontsize=7,fontweight="normal",color="#444"); lab(ax,"A",-.25)
# B: candidates heatmap
ax=fig.add_subplot(gsp[:,1]); X=lg.loc[cand.gene,ordl]; Zs=X.sub(X.mean(axis=1),axis=0).div(X.std(axis=1),axis=0)
im=ax.imshow(Zs.values,cmap="RdBu_r",vmin=-2.2,vmax=2.2,aspect="auto")
ax.set_yticks(range(len(cand))); ax.set_yticklabels([f"{g} (r={r:.2f})" for g,r in zip(cand.gene,cand.hipsci_rna_protein_corr)],fontsize=5.2,style="italic")
ax.set_xticks(range(14)); ax.set_xticklabels([l.replace("HEL","") for l in ordl],rotation=90,fontsize=5.8)
for j,l in enumerate(ordl): ax.add_patch(plt.Rectangle((j-.5,-1.5),1,.8,color=C[la.loc[l,"group"]],clip_on=False))
ax.set_ylim(len(cand)-.5,-1.6); ax.spines[:].set_visible(False); ax.tick_params(length=0)
cb=fig.colorbar(im,ax=ax,fraction=.04,pad=.02); cb.set_label("z-score",fontsize=6.5)
ax.set_title(f"{len(cand)} genes: nominal p < 0.05 and HipSci RNA–protein r ≥ 0.3\n(none at FDR < 0.1; no enriched GO term)",fontsize=7,fontweight="normal",color="#444",pad=12); lab(ax,"B",-.38)
# C: what the iPSC proteome captures / misses
ax=fig.add_subplot(gsp[1,0])
pick=[("expressed, protein detected in HipSci (8,726)",["Metabolism of RNA","Translation","Spliceosome","mitochondrial matrix","Cell Cycle"]),
      ("expressed, protein not detected (4,514)",["G protein-coupled receptor activity","monoatomic ion channel activity","extracellular matrix","Cytokine-cytokine receptor interaction","transmembrane signaling receptor activity"])]
rows=[]
for s,terms in pick:
    for t in terms:
        r=go[(go.gene_set==s)&(go.term_name==t)].sort_values("p_adj_gSCS").iloc[0]
        rows.append((s,t,r.fold_enrichment,r.p_adj_gSCS))
rows=rows[::-1]; yy=np.arange(len(rows))
cols=["#2a78d6" if "not detected" not in s else "#9aa0a6" for s,_,_,_ in rows]
ax.barh(yy,[f for _,_,f,_ in rows],color=cols,height=.65)
ax.axvline(1,color=MUTED,lw=.7,ls="--")
ax.set_yticks(yy); ax.set_yticklabels([t[0].upper()+t[1:] if t[0].islower() else t for _,t,_,_ in rows],fontsize=5.8)
ax.set_xlabel("Fold enrichment vs expressed protein-coding genes"); ax.set_xlim(0,2.6)
for yi,(_,_,f,p) in zip(yy,rows): ax.text(f+.04,yi,"p<1e-10" if p<1e-10 else f"p={p:.0e}",va="center",fontsize=5.2,color="#444")
ax.legend(handles=[Patch(color="#2a78d6",label="Protein detected"),Patch(color="#9aa0a6",label="Protein not detected")],
          frameon=False,fontsize=5.8,loc="upper center",bbox_to_anchor=(0.4,-0.13),ncol=2)
ax.set_title("GO / pathway enrichment (g:Profiler, g:SCS)",fontsize=7,fontweight="normal",color="#444"); lab(ax,"C",-.9)
save(fig,"FigureS3_protein_crossref")
print("done")
