# R scripts for the main analyses

Run from this folder (RStudio: Session > Set Working Directory > To Source File Location). Input: `../../data/analysis_dataset_110_patients.csv`.
Packages: survival, glmnet, ggplot2, parallel, pec, rms, timeROC (optional cross-check). R 4.6.0 (`../../results/revision_analyses_output/sessionInfo.txt`).

| Script | Content | Manuscript item |
|---|---|---|
| `10_tace_revision_analyses.R` | **Main script of the revision analyses** (sections 0-10): reproduction of Tables 3 and 5, apparent comparison, bootstrap optimism with selection repeated, 5 x 10-fold CV with LASSO / elastic net, 70/30 splits, feature-family ablation, extended Cox model, decision curves, size adjustment, Table 1. Full run: 15-40 min (`TACE_TEST=1` for a 5-minute test; `TACE_N_CORES` sets the cores). Writes to `../../results/revision_analyses_output/` | Tables 3-5, Figs 4, 6, 7, Supplementary Tables S2-S10 |
| `00_common.R`, `01_univariate_fdr_stepwise.R` | Univariate Cox with Benjamini-Hochberg correction and backward elimination of the full cohort | Table 2, Supplementary Table S3 |
| `02_models_ph_extended_cox.R` | Tumour-burden adjustment, Schoenfeld tests, extended Cox model, time-specific and interval HRs | Supplementary Tables S4, S5 |
| `05_brier_and_calibration_slope.R` | IPCW Brier scores; bootstrap calibration slope | Calibration section |
| `11_make_fig5_roc.R` | Fig 5 (time-dependent ROC curves, fixed three-feature model vs HAP) | Fig 5, Table 4 AUCs |
| `../00_run_checks.R` | timeROC AUCs and bootstrap-corrected C-index (`rms::validate`) | Table 4 |

Selection rule note: in every selection-aware analysis the pipeline is univariate Cox + Benjamini-Hochberg (q < 0.05) -> backward elimination (P < 0.05, HAP and BCLC forced in); if more than 25 features pass FDR screening, the 25 with the smallest P values enter backward elimination (numerical safeguard).
