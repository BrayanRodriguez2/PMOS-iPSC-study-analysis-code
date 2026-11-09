import openpyxl, numpy as np
from scipy import stats, optimize
wb = openpyxl.load_workbook("Additional_file_1.xlsx", read_only=True)
rows = list(wb["S4_DE"].iter_rows(values_only=True))[1:]
fc = np.array([r[1] for r in rows], float); p = np.array([r[2] for r in rows], float)
m = (p < 0.95) & (p > 0)
t = stats.t.isf(p[m] / 2, 12); sd = np.abs(fc[m]) / t / np.sqrt(2 / 7)
SD = np.median(sd)
def mdd(n, alpha, power=0.8):
    df = 2 * n - 2; tc = stats.t.isf(alpha / 2, df)
    return optimize.brentq(lambda d: stats.nct.sf(tc, df, d / SD * np.sqrt(n / 2)) - power, 1e-3, 20)
for n in (3, 4, 7, 10, 15, 20):
    print(n, round(mdd(n, 0.05), 2), round(mdd(n, 0.001), 2))
