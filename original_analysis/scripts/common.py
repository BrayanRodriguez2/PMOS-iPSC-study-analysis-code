"""Shared helpers for the Python steps (imported by 05-09)."""
import re
import numpy as np, pandas as pd
from scipy import stats
 
 
# The RNA-seq files label the control clone of donor 114 as HEL362_3; the correct
# clone is HEL362.1 (UHi034-A), confirmed 2 Nov 2026.
LINE_FIX = {"HEL362.3": "HEL362.1"}
 
 
def hel_name(x: str) -> str:
    """'HEL358_1_S8_L001_R1_001' / 'HEL 358.1' / 'HEL358.1' -> 'HEL358.1'"""
    m = re.search(r"HEL\s*_?(\d+)[._](\d+)", str(x))
    name = f"HEL{m.group(1)}.{m.group(2)}" if m else str(x)
    return LINE_FIX.get(name, name)
 
 
def read_key(path):
    k = pd.read_csv(path)
    k["array_iid"] = k["array_iid"].astype(str)
    k["hel_line"] = k["hel_line"].fillna("").map(lambda s: hel_name(s) if s else "")
    return k
 
 
def read_counts(path):
    c = pd.read_csv(path, sep="\t", index_col=0)
    c.columns = [hel_name(x) for x in c.columns]
    return c
 
 
def ols(y, X):
    """OLS with intercept. Returns (beta, se, t, p) for every column of X."""
    X = np.column_stack([np.ones(len(y)), X])
    n, k = X.shape
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ b
    df = n - k
    s2 = resid @ resid / df
    cov = s2 * np.linalg.pinv(X.T @ X)
    se = np.sqrt(np.diag(cov))
    t = b / se
    p = 2 * stats.t.sf(np.abs(t), df)
    return b[1:], se[1:], t[1:], p[1:]
 
 
def bh(p):
    p = np.asarray(p, float); n = len(p)
    o = np.argsort(p); r = np.empty(n); r[o] = np.arange(1, n + 1)
    q = p * n / r
    qs = np.minimum.accumulate(q[o][::-1])[::-1]
    out = np.empty(n); out[o] = np.minimum(qs, 1)
    return out
 
 
def residualize(y, covars):
    """Residuals of y after regressing on covariates (keeps the mean)."""
    X = np.column_stack([np.ones(len(y)), covars])
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    return y - X @ b + y.mean()
