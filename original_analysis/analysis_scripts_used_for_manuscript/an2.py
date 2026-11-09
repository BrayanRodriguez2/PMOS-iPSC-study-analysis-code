import numpy as np, pandas as pd
from scipy import stats
from itertools import combinations
lg=pd.read_csv("/home/claude/ms/tables/log2_expressed.tsv.gz",sep="\t",index_col=0)
key=pd.read_csv("/home/claude/pheno/sample_key_14_ipsc_donors.csv").set_index("hel_line").loc[lg.columns]
cnt_x=None
R=np.apply_along_axis(stats.rankdata,1,lg.values); R=R-R.mean(1,keepdims=True); R/=np.sqrt((R**2).sum(1,keepdims=True))
n=R.shape[1]
def nsig(x):
    rx=stats.rankdata(x); rx=rx-rx.mean(); rx/=np.sqrt((rx**2).sum())
    r=R@rx; t=r*np.sqrt((n-2)/(1-r**2+1e-12)); p=2*stats.t.sf(np.abs(t),n-2)
    return (p<0.05).sum()
rng=np.random.default_rng(42)
xist=np.log2(pd.read_csv("/mnt/user-data/attachments/522de043-d88c-536e-9f44-66e7182b8053/336ebae1-Normalized_counts.tabular",sep="\t",index_col=0).filter(like="HEL").loc["XIST"].values+1)
traits={"PMOS status":(key.group=="PMOS").astype(float).values,"AMH":key.amh.values,"Testosterone":key.testosterone.values,"FAI":key.fai.values,
        "FG score":key.fg_score.values,"BMI":key.bmi.values,"HOMA-IR":key.homa_ir.values,"Age":key.age.values}
out=[]
for k,x in traits.items():
    obs=nsig(x)
    if k=="PMOS status":
        null=[]
        for c in combinations(range(n),7):
            z=np.zeros(n); z[list(c)]=1; null.append(nsig(z))
        null=np.array(null)
    else:
        null=np.array([nsig(rng.permutation(x)) for _ in range(2000)])
    out.append(dict(trait=k,genes_p05=int(obs),null_median=int(np.median(null)),null_95pct=int(np.percentile(null,95)),perm_p=(np.sum(null>=obs)+1)/(len(null)+1),n_perm=len(null)))
    print(out[-1],flush=True)
pd.DataFrame(out).to_csv("/home/claude/ms/tables/S5b_global_permutation.tsv",sep="\t",index=False)
