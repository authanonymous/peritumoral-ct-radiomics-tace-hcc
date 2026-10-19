# =============================================================================
# 01_univariate_fdr_stepwise.R
# Univariate Cox screening of the 624 features with Benjamini-Hochberg FDR, then multivariable Cox with HAP score and
# BCLC stage mandatory and backward elimination (P < 0.05) of the radiomic candidates -> fixed three-feature model.
# Reproduces: Table 2, Supplementary Table S3.   Outputs: ../../results/table_S3_univariate_cox_all.csv, table_2_final_model.csv
# Expected (manuscript): 5 features with q < 0.05; final HR (95% CI): MaxProb 1.316 (1.097-1.579), MCC 1.681 (1.276-2.216),
#   Imc1 1.586 (1.187-2.119), HAP 1.602 (1.212-2.116), BCLC 1.601 (0.970-2.645).
# =============================================================================
source("00_common.R")
d <- load_data(); feats <- feature_names(d)
cat("Patients:", nrow(d), " Deaths:", sum(d$event), " Features:", length(feats), "\n\n")

scr <- uni_screen(d, feats)
scr <- scr[order(scr$p_value), ]
write.csv(scr, file.path(OUT, "table_S3_univariate_cox_all.csv"), row.names = FALSE)
cat("Features with FDR q < 0.05:", sum(scr$FDR_q < 0.05), " (not retained:", sum(scr$FDR_q >= 0.05), ")\n")
print(scr[scr$FDR_q < 0.05, c("Feature", "HR", "CI_low", "CI_high", "p_value", "FDR_q")], digits = 4, row.names = FALSE)

sp <- select_pipeline(d, feats)
cat("\nSelected after backward elimination:\n"); print(sp$features)
sm <- summary(sp$fit)$coefficients; ci <- exp(confint(sp$fit))
tab <- data.frame(Variable = gsub("`", "", rownames(sm), fixed = TRUE), HR = sm[, "exp(coef)"],
                  CI_low = ci[, 1], CI_high = ci[, 2], P = sm[, "Pr(>|z|)"], row.names = NULL)
write.csv(tab, file.path(OUT, "table_2_final_model.csv"), row.names = FALSE)
cat("\nTable 2 (final multivariable model):\n"); print(tab, digits = 4, row.names = FALSE)

cat("\nCorrelation among the selected features (Pearson):\n")
print(round(cor(d[, sp$features]), 2))
saveRDS(sp$features, file.path(OUT, "selected_features_full_cohort.rds"))
stopifnot(setequal(sp$features, FIXED3))   # the three features used as the "fixed three-feature model"
