# ------------------------------------------------------------------
# 00_run_checks.R  -  time-dependent AUC (timeROC) and bootstrap-corrected
# C-index (rms::validate) for the fixed three-feature model, HAP+BCLC and HAP.
# Run from the folder "code/" (RStudio: Session > Set Working Directory > To Source File Location).
# Packages: survival, timeROC, rms   install.packages(c("survival","timeROC","rms"))
# Output reported in results/run_checks_output.txt
# ------------------------------------------------------------------
library(survival); library(timeROC); library(rms)

DATA <- "../data/analysis_dataset_110_patients.csv"
if (!file.exists(DATA)) stop("Run this script from the code/ folder; missing: ", DATA)

df <- read.csv(DATA, check.names = FALSE)
df$os <- df$survival_time / 30.44                     # days -> months
df$z1 <- df[["ratio_0-3_vs_5-10_glcm_MaximumProbability"]]
df$z2 <- df[["ratio_intra_vs_0-3_glcm_MCC"]]
df$z3 <- df[["ratio_3-5_vs_5-10_glcm_Imc1"]]
cat("Patients:", nrow(df), " Events:", sum(df$death), "\n\n")

forms <- list(
  "Fixed three-feature model" = Surv(os, death) ~ z1 + z2 + z3 + hap_score + bclc,
  "HAP+BCLC"                  = Surv(os, death) ~ hap_score + bclc,
  "HAP alone"                 = Surv(os, death) ~ hap_score)

cat("=== 1) timeROC AUC at 1, 2, 3 years ===\n")
for (nm in names(forms)) {
  fit <- coxph(forms[[nm]], data = df)
  lp  <- predict(fit, type = "lp")
  r   <- timeROC(T = df$os, delta = df$death, marker = lp, cause = 1,
                 weighting = "marginal", times = c(12, 24, 36), iid = FALSE)
  cat(sprintf("%-26s apparent C = %.3f | AUC 1/2/3 y = %.3f / %.3f / %.3f\n",
              nm, concordance(fit)$concordance, r$AUC[1], r$AUC[2], r$AUC[3]))
}

cat("\n=== 2) Bootstrap-corrected C-index (rms::validate, B = 500) ===\n")
for (nm in names(forms)) {
  f <- cph(forms[[nm]], data = df, x = TRUE, y = TRUE)
  set.seed(42)
  v <- validate(f, method = "boot", B = 500)
  cat(sprintf("%-26s apparent C = %.3f | optimism = %.3f | corrected C = %.3f\n", nm,
      0.5 + v["Dxy", "index.orig"] / 2,
      (v["Dxy", "optimism"]) / 2,
      0.5 + v["Dxy", "index.corrected"] / 2))
}
