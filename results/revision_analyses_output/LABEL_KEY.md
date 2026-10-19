# Label key for the output tables in this folder

The script `code/R_analysis/10_tace_revision_analyses.R` uses short internal labels. They correspond to the names used in the manuscript and supplement as follows.

| Label in output files | Name in the manuscript |
|---|---|
| M1_HAP | HAP alone |
| M2_HAP_BCLC | HAP+BCLC |
| M3_fixed_3feat (also "M3 (published)", `M3_fixed_3feat_LEAKY`) | Fixed three-feature model (features chosen once on the full cohort, so its out-of-sample performance is leakage-prone) |
| M4_stepwise, M4_selection_on_train, M4_nested_pipe | Stepwise radiomic pipeline with selection repeated in every training set or resample |
| LASSO_min / LASSO_1se | LASSO Cox, minimum-deviance / 1-SE penalty |
| EN0.5_min / EN0.5_1se | Elastic net (alpha = 0.5) Cox, minimum-deviance / 1-SE penalty |
| `published_feature` (TRUE/FALSE) | Feature belongs to the fixed three-feature model |
| ratio_0.3_vs_5.10_glcm_MaximumProbability | GLCM Maximum Probability ratio (0-3 / 5-10 mm) |
| ratio_intra_vs_0.3_glcm_MCC | GLCM MCC ratio (intratumoral / 0-3 mm) |
| ratio_3.5_vs_5.10_glcm_Imc1 | GLCM Imc1 ratio (3-5 / 5-10 mm) |
| 01_... to 10_... | Section numbers of the script (01 Table 3/5 reproduction, 02 apparent models, 03 bootstrap, 04 cross-validation, 05 splits, 06 ablation, 07 extended Cox, 08 decision curves, 09 size adjustment, 10 Table 1) |

Run with R 4.6.0 (see `sessionInfo.txt`); seed 20260920.
