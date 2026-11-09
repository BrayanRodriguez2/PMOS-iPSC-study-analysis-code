import sys; sys.path.insert(0,"/home/claude/pcos_ipsc_genetics/scripts")
import numpy as np, pandas as pd
from scipy import stats
from common import read_counts, ols, bh
OUT="/home/claude/ms"
cnt=read_counts("/mnt/user-data/attachments/522de043-d88c-536e-9f44-66e7182b8053/336ebae1-Normalized_counts.tabular")
key=pd.read_csv("/home/claude/pheno/sample_key_14_ipsc_donors.csv",dtype={"array_iid":str}).set_index("hel_line")
ph=pd.read_csv("/home/claude/pheno/donor_phenotypes.tsv",sep="\t")
lines=list(key.sort_values(["group","bmi"]).index); cnt=cnt[lines]
expd=cnt[(cnt>=10).sum(1)>=7]            # expressed: >=10 in at least half the lines
lg=np.log2(expd+1)
grp=(key.loc[lines,"group"]=="PMOS").astype(float).values
xist=np.log2(cnt.loc["XIST",lines]+1).values
print("expressed genes",len(lg))
# DE models
res={}
for lab,X in [("group",grp[:,None]),("group+XIST",np.column_stack([grp,xist]))]:
    B=[];P=[]
    for g,y in zip(lg.index,lg.values):
        b,se,t,p=ols(y,X); B.append(b[0]);P.append(p[0])
    res[lab]=pd.DataFrame({"gene":lg.index,"log2FC_PMOS_vs_control":B,"p":P}); res[lab]["q_BH"]=bh(P)
    r=res[lab]; print(lab,"p<0.05:",(r.p<0.05).sum(),"p<0.001:",(r.p<0.001).sum(),"q<0.1:",(r.q_BH<0.1).sum(),"q<0.05:",(r.q_BH<0.05).sum())
de=res["group"].merge(res["group+XIST"],on="gene",suffixes=("","_adjXIST"))
de["mean_norm_count"]=expd.mean(1).values
de.sort_values("p").to_csv(f"{OUT}/tables/S4_DE_PMOS_vs_control.tsv",sep="\t",index=False)
# phenotype associations
phen={"AMH":"amh","Testosterone":"testosterone","FAI":"fai","FG score":"fg_score","BMI":"bmi","HOMA-IR":"homa_ir","Age":"age"}
assoc=[]; summ=[]
for name,c in phen.items():
    x=key.loc[lines,c].values.astype(float)
    r=np.array([stats.spearmanr(x,y)[0] for y in lg.values]); p=np.array([stats.spearmanr(x,y)[1] for y in lg.values])
    q=bh(p); summ.append(dict(trait=name,genes_p05=int((p<0.05).sum()),genes_q10=int((q<0.1).sum()),genes_q05=int((q<0.05).sum()),expected_p05=round(0.05*len(p))))
    assoc.append(pd.DataFrame({"gene":lg.index,"trait":name,"rho":r,"p":p,"q_BH":q}))
for name,x in [("XIST (log2)",xist),("PMOS status",grp)]:
    r=np.array([stats.spearmanr(x,y)[0] for y in lg.values]); p=np.array([stats.spearmanr(x,y)[1] for y in lg.values]); q=bh(p)
    summ.append(dict(trait=name,genes_p05=int((p<0.05).sum()),genes_q10=int((q<0.1).sum()),genes_q05=int((q<0.05).sum()),expected_p05=round(0.05*len(p))))
    assoc.append(pd.DataFrame({"gene":lg.index,"trait":name,"rho":r,"p":p,"q_BH":q}))
A=pd.concat(assoc); S=pd.DataFrame(summ); print(S.to_string(index=False))
A[A.p<0.01].sort_values(["trait","p"]).to_csv(f"{OUT}/tables/S5_gene_trait_associations_p01.tsv",sep="\t",index=False)
S.to_csv(f"{OUT}/tables/S5a_trait_association_summary.tsv",sep="\t",index=False)
# PCA
v=lg.var(1).sort_values(ascending=False).index[:1000]
Z=lg.loc[v].sub(lg.loc[v].mean(1),axis=0).T.values
U,Sv,Vt=np.linalg.svd(Z,full_matrices=False); pcs=U*Sv; ve=Sv**2/np.sum(Sv**2)
pcdf=pd.DataFrame(pcs[:,:5],index=lines,columns=[f"PC{i}" for i in range(1,6)])
pcdf.to_csv(f"{OUT}/tables/pcs.tsv",sep="\t"); print("PC var%",np.round(ve[:5]*100,1))
cov={"PMOS":grp,"XIST":xist,**{k:key.loc[lines,c].values.astype(float) for k,c in phen.items()}}
pcc=pd.DataFrame({k:[stats.spearmanr(val,pcdf[f"PC{i}"])[0] for i in range(1,6)] for k,val in cov.items()},index=[f"PC{i} ({ve[i-1]*100:.0f}%)" for i in range(1,6)])
pcp=pd.DataFrame({k:[stats.spearmanr(val,pcdf[f"PC{i}"])[1] for i in range(1,6)] for k,val in cov.items()},index=pcc.index)
pcc.to_csv(f"{OUT}/tables/pc_trait_rho.tsv",sep="\t"); pcp.to_csv(f"{OUT}/tables/pc_trait_p.tsv",sep="\t")
print(pcc.round(2).to_string()); print(pcp.round(3).to_string())
# variance explained per factor (univariate R2)
R2={}
for k,val in cov.items():
    xv=(val-val.mean()); 
    yy=lg.values-lg.values.mean(1,keepdims=True)
    r=(yy@xv)/np.sqrt((yy**2).sum(1)*(xv**2).sum()); R2[k]=r**2
R2=pd.DataFrame(R2,index=lg.index); R2.to_csv(f"{OUT}/tables/univariate_R2.tsv.gz",sep="\t")
print(R2.median().round(3).to_string())
# XIST correlated genes
xr=A[A.trait=="XIST (log2)"].sort_values("p").head(30); print(xr[["gene","rho","p","q_BH"]].to_string(index=False))
lg.to_csv(f"{OUT}/tables/log2_expressed.tsv.gz",sep="\t")
