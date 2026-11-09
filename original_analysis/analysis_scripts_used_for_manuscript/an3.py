import numpy as np, pandas as pd
from scipy import stats
cnt=pd.read_csv("/mnt/user-data/attachments/522de043-d88c-536e-9f44-66e7182b8053/336ebae1-Normalized_counts.tabular",sep="\t",index_col=0)
cnt.columns=[ "HEL"+c.split("_")[0][3:]+"."+c.split("_")[1] for c in cnt.columns]
lg=pd.read_csv("/home/claude/ms/tables/log2_expressed.tsv.gz",sep="\t",index_col=0)
key=pd.read_csv("/home/claude/pheno/sample_key_14_ipsc_donors.csv").set_index("hel_line").loc[lg.columns]
st=pd.read_csv("/home/claude/geno/leads29.tsv",sep="\t")
gmap={"HHAT":["HHAT"],"THADA":["THADA"],"MSH6":["MSH6"],"FSHR":["FSHR"],"INHBB":["INHBB"],"TEX41":["TEX41"],"ERBB4":["ERBB4"],
"TOPAZ1":["TOPAZ1"],"MYO10":["MYO10"],"IRF / RAD50":["IRF1","RAD50"],"ID4":["ID4"],"CYP3 complex":["CYP3A7","CYP3A5"],"NEIL2 / GATA4":["NEIL2","GATA4"],
"FANCC":["FANCC"],"DENND1A":["DENND1A"],"FSHB":["FSHB"],"PRSS23":["PRSS23"],"YAP1":["YAP1"],"ZBTB16":["ZBTB16"],"ERBB3":["ERBB3"],
"HELB / HMGA2":["HELB","HMGA2"],"KRR1":["KRR1"],"TOX3":["TOX3"],"FTO":["FTO"],"MAF":["MAF"],"SHBG":["SHBG"],"AMH":["AMH"],"MAPRE1":["MAPRE1"],"CHEK2":["CHEK2"]}
rows=[]
for _,r in st.iterrows():
    for g in gmap[r.consensus]:
        m=cnt.loc[g,lg.columns] if g in cnt.index else None
        rows.append(dict(locus_rsid=r.rsid,consensus=r.consensus,gene=g,mean_norm_count=None if m is None else round(m.mean(),1),
                         expressed= g in lg.index))
G=pd.DataFrame(rows)
# percentile of expression among all expressed genes
allmean=cnt.loc[lg.index,lg.columns].mean(1)
G["expression_percentile"]=[round(stats.percentileofscore(allmean,v),0) if e else np.nan for v,e in zip(G.mean_norm_count,G.expressed)]
de=pd.read_csv("/home/claude/ms/tables/S4_DE_PMOS_vs_control.tsv",sep="\t").set_index("gene")
G=G.join(de[["log2FC_PMOS_vs_control","p","q_BH","log2FC_PMOS_vs_control_adjXIST","p_adjXIST"]],on="gene")
for t,c in [("AMH","amh"),("FAI","fai"),("BMI","bmi")]:
    G[f"rho_{t}"]=[stats.spearmanr(key[c],lg.loc[g])[0] if g in lg.index else np.nan for g in G.gene]
    G[f"p_{t}"]=[stats.spearmanr(key[c],lg.loc[g])[1] if g in lg.index else np.nan for g in G.gene]
G.to_csv("/home/claude/ms/tables/S7_PCOS_GWAS_genes_expression.tsv",sep="\t",index=False)
print(G[["gene","mean_norm_count","expressed","expression_percentile","log2FC_PMOS_vs_control","p","rho_AMH","p_AMH","rho_FAI","p_FAI"]].round(3).to_string(index=False))
print("expressed:",G.expressed.sum(),"of",len(G))
# gene-set score tests (expression-matched permutation)
gs=[g for g in G.gene[G.expressed]]
Z=lg.sub(lg.mean(1),axis=0).div(lg.std(1),axis=0)
score=Z.loc[gs].mean()
bins=pd.qcut(lg.mean(1),20,labels=False); need=bins.loc[gs].value_counts()
rng=np.random.default_rng(7); pos={b:np.where(bins.values==b)[0] for b in need.index}
traits={"PMOS status":(key.group=="PMOS").astype(float),"AMH":key.amh,"FAI":key.fai,"Testosterone":key.testosterone,"BMI":key.bmi,"HOMA-IR":key.homa_ir}
nulls={k:[] for k in traits}
Zv=Z.values
for i in range(10000):
    pick=np.concatenate([rng.choice(pos[b],n,replace=False) for b,n in need.items()])
    s=Zv[pick].mean(0)
    for k,x in traits.items(): nulls[k].append(stats.spearmanr(x,s)[0])
out=[]
for k,x in traits.items():
    r=stats.spearmanr(x,score)[0]; nl=np.array(nulls[k])
    out.append(dict(test=f"PCOS GWAS gene-set score vs {k}",n_genes=len(gs),rho=round(r,3),p_perm_matched=(np.sum(np.abs(nl)>=abs(r))+1)/(len(nl)+1)))
O=pd.DataFrame(out); print(O.to_string(index=False)); O.to_csv("/home/claude/ms/tables/S7b_geneset_tests.tsv",sep="\t",index=False)
pd.DataFrame({"score":score,"group":key.group}).to_csv("/home/claude/ms/tables/geneset_score.tsv",sep="\t")
