# Results folder

- `revision_analyses_output/`: complete output of `code/R_analysis/10_tace_revision_analyses.R` (R 4.6.0, seed 20260920): bootstrap with selection repeated, repeated cross-validation, 100 splits, feature-family ablation, extended Cox, decision curves, size adjustment, Table 1; `sessionInfo.txt`; `LABEL_KEY.md` maps output labels to the manuscript names.
- `S1_by_exclusion_group.csv`: Supplementary Table S1 (analysed patients vs. each exclusion group).
- `preprocessing_step_counts.csv`: feature counts at each preprocessing step.
- `run_checks_output.txt`: output of `code/00_run_checks.R` (timeROC AUC and rms bootstrap).
Running the scripts in `code/R_analysis/` 01, 02 and 05 writes further tables here (`table_2_final_model.csv`, `table_S3_univariate_cox_all.csv`, `table_S4_adjusted_models.csv`, `table_S5_time_varying.csv`, `table_brier_ipcw.csv`, `table_calibration_slope_bootstrap.csv`).
