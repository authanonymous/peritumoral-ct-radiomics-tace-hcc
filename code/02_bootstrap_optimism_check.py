"""
Independent Python cross-check of the conventional (Harrell) optimism bootstrap.
Shows (1) the corrected routine: rows resampled together (covariates AND outcome) and
          (2) a replica of the original faulty routine in which the outcome was NOT resampled,
which reproduces the earlier reported values (HAP alone: optimism ~0.141, corrected C ~0.51;
fixed three-feature model: ~0.64).  Deposit your own R scripts as the primary code; this file is a cross-check.

Run from code/ . Input: ../data/analysis_dataset_110_patients.csv
pip install pandas numpy lifelines
"""
import numpy as np, pandas as pd, warnings
warnings.filterwarnings("ignore")
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index

df = pd.read_csv("../data/analysis_dataset_110_patients.csv")
d = pd.DataFrame({"z1": df["ratio_0-3_vs_5-10_glcm_MaximumProbability"], "z2": df["ratio_intra_vs_0-3_glcm_MCC"],
                  "z3": df["ratio_3-5_vs_5-10_glcm_Imc1"], "hap": df["hap_score"].astype(float),
                  "bclc": df["bclc"].astype(float), "t": df["survival_time"].astype(float) / 30.44, "e": df["death"].astype(int)})
models = {"Fixed three-feature model": ["z1", "z2", "z3", "hap", "bclc"], "HAP+BCLC": ["hap", "bclc"], "HAP": ["hap"]}
fit = lambda data, cols: CoxPHFitter().fit(data[cols + ["t", "e"]], "t", "e")
C = lambda m, data, cols: concordance_index(data.t, -m.predict_log_partial_hazard(data[cols]).values, data.e)
rng = np.random.default_rng(42)
for name, cols in models.items():
    m = fit(d, cols); app = C(m, d, cols)
    good = []
    for _ in range(500):
        b = d.iloc[rng.integers(0, len(d), len(d))].reset_index(drop=True)      # outcome resampled WITH covariates
        try: mb = fit(b, cols); good.append(C(mb, b, cols) - C(mb, d, cols))
        except Exception: pass
    bad = []
    for _ in range(1000):
        idx = rng.integers(0, len(d), len(d))
        b = d[cols].iloc[idx].reset_index(drop=True); b["t"] = d.t.values; b["e"] = d.e.values   # outcome NOT resampled (the bug)
        try:
            mb = fit(b, cols); o = C(mb, b, cols) - concordance_index(d.t, mb.predict_log_partial_hazard(d[cols]).values, d.e)
            if o > -0.1: bad.append(o)
        except Exception: pass
    print(f"{name:28s} apparent {app:.3f} | CORRECT optimism {np.mean(good):.3f} -> corrected {app-np.mean(good):.3f} (n={len(good)})"
          f" | FAULTY replica optimism {np.mean(bad):.3f} -> {app-np.mean(bad):.3f} (valid {len(bad)}/1000)")
