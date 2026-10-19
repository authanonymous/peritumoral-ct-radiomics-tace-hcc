"""
Refits the extended Cox model (GLCM Maximum Probability ratio x log(t+1), months) with lifelines, checks it against
Table 3B / Supplementary Table S5 (interaction HR 1.64, 1.13-2.39; main HR 0.61, 0.31-1.22) and redraws Fig 4 with
two at-risk thresholds (<25 and <10 patients at risk, computed from the data).
Run from code/ . Input: ../data/analysis_dataset_110_patients.csv   pip install pandas numpy lifelines matplotlib
"""
import numpy as np, pandas as pd, warnings
warnings.filterwarnings("ignore")
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from lifelines import CoxTimeVaryingFitter
df = pd.read_csv("../data/analysis_dataset_110_patients.csv")
d = pd.DataFrame({"id": range(len(df)), "z1": df["ratio_0-3_vs_5-10_glcm_MaximumProbability"], "z2": df["ratio_intra_vs_0-3_glcm_MCC"],
                  "z3": df["ratio_3-5_vs_5-10_glcm_Imc1"], "hap": df["hap_score"].astype(float), "bclc": df["bclc"].astype(float),
                  "t": df["survival_time"].astype(float) / 30.44, "e": df["death"].astype(int)})
ev = np.sort(d.t[d.e == 1].unique()); rows = []
for _, r in d.iterrows():                       # split follow-up at every event time (equivalent to R tt())
    prev = 0.0
    for te in ev:
        if te > r.t: break
        rows.append((r.id, prev, te, r.z1, r.z2, r.z3, r.hap, r.bclc, r.z1 * np.log(te + 1), int(r.e == 1 and te == r.t))); prev = te
    if not (r.e == 1 and prev == r.t) and r.t > prev:
        rows.append((r.id, prev, r.t, r.z1, r.z2, r.z3, r.hap, r.bclc, r.z1 * np.log(r.t + 1), 0))
L = pd.DataFrame(rows, columns=["id", "start", "stop", "z1", "z2", "z3", "hap", "bclc", "z1xlog", "event"])
m = CoxTimeVaryingFitter().fit(L, id_col="id", event_col="event", start_col="start", stop_col="stop")
print(np.exp(m.summary[["coef", "coef lower 95%", "coef upper 95%"]]).round(3))
names = list(m.params_.index); V = np.asarray(m.variance_matrix_); i, j = names.index("z1"), names.index("z1xlog")
b1, b2 = float(m.params_["z1"]), float(m.params_["z1xlog"])
t = np.linspace(1, 96, 400); l = np.log(t + 1); lh = b1 + b2 * l; se = np.sqrt(V[i, i] + l * l * V[j, j] + 2 * l * V[i, j])
hr, lo, hi = np.exp(lh), np.exp(lh - 1.96 * se), np.exp(lh + 1.96 * se)
grid = np.arange(0, 100, 0.25); atrisk = np.array([(d.t >= g).sum() for g in grid])
T25, T10 = grid[np.argmax(atrisk < 25)], grid[np.argmax(atrisk < 10)]; print("thresholds (months):", T25, T10)
fig, ax = plt.subplots(figsize=(10.5, 6.75), dpi=200)
def band(a, b, c, al): k = (t >= a) & (t <= b); ax.fill_between(t[k], lo[k], hi[k], color=c, alpha=al, lw=0)
band(0, T25, "#7F9DB5", .45); band(T25, T10, "#B0B0B0", .55); band(T10, 100, "#6E6E6E", .55)
for a, b, ls in ((0, T25, ":"), (T25, T10, "--"), (T10, 100, "-")):
    k = (t >= a) & (t <= b); ax.plot(t[k], hr[k], "k", lw=1.6, ls=ls)
ax.axhline(1, color="#555", ls=(0, (5, 4))); ax.set_yscale("log"); ax.set_yticks([1, 3, 10]); ax.set_yticklabels(["1", "3", "10"])
ax.set_ylim(.6, 25); ax.set_xlim(-3, 100); ax.set_xlabel("Months since TACE"); ax.set_ylabel("HR per SD, GLCM MaxProb ratio (log scale)")
ax.set_title("Time-varying HR (extended Cox model)", loc="left"); ax.grid(color="#E6E6E6")
ax.legend(handles=[Patch(fc="#7F9DB5", alpha=.6, label="95% CI, \u2265 25 patients at risk"), Patch(fc="#B0B0B0", alpha=.8, label="< 25 patients at risk"),
                   Patch(fc="#6E6E6E", alpha=.8, label="< 10 patients at risk")], loc="upper left", fontsize=7, frameon=False)
fig.tight_layout(); fig.savefig("../results/Fig4_revised.png"); print("saved ../results/Fig4_revised.png")
