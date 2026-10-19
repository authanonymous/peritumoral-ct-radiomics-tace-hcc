# =============================================================================
# 05_brier_and_calibration_slope.R
# (a) IPCW Brier scores at 1, 2 and 3 years (apparent; Kaplan-Meier censoring model) for HAP, HAP+BCLC and the fixed three-feature model
# (b) bootstrap-estimated calibration slope of the linear predictor (rms::validate, 500 resamples; fixed feature set)
# Reproduces: "Calibration and Prediction Error" (Fig 6).   Outputs: ../../results/table_brier_ipcw.csv, table_calibration_slope_bootstrap.csv
# Expected (manuscript): Brier 0.144 / 0.201 / 0.197 (fixed model), 0.162 / 0.213 / 0.219 (HAP+BCLC), 0.162 / 0.213 / 0.224 (HAP);
#   bootstrap calibration slope 1.04 (HAP), 0.96 (HAP+BCLC), 0.78 (fixed three-feature model).
# Packages: survival, pec, rms
# =============================================================================
source("00_common.R")
for (p in c("pec", "rms")) if (!requireNamespace(p, quietly = TRUE)) stop("Please install.packages('", p, "')")
suppressPackageStartupMessages({ library(pec); library(rms) })
d <- load_data()
d$z1 <- d[["ratio_0-3_vs_5-10_glcm_MaximumProbability"]]; d$z2 <- d[["ratio_intra_vs_0-3_glcm_MCC"]]; d$z3 <- d[["ratio_3-5_vs_5-10_glcm_Imc1"]]
TIMES <- c(12, 24, 36)

# ---- (a) IPCW Brier score
m_hap <- coxph(Surv(os_months, event) ~ hap,                          data = d, x = TRUE, y = TRUE)
m_hb  <- coxph(Surv(os_months, event) ~ hap + bclc,                   data = d, x = TRUE, y = TRUE)
m_fix <- coxph(Surv(os_months, event) ~ z1 + z2 + z3 + hap + bclc,    data = d, x = TRUE, y = TRUE)
pe <- pec(list(`Fixed three-feature model` = m_fix, `HAP+BCLC` = m_hb, HAP = m_hap), formula = Surv(os_months, event) ~ 1, data = d,
          times = TIMES, cens.model = "marginal", splitMethod = "none", reference = TRUE, verbose = FALSE)
brier <- sapply(names(pe$AppErr), function(n) sapply(TIMES, function(tm) pe$AppErr[[n]][which.min(abs(pe$time - tm))]))
rownames(brier) <- paste0(TIMES / 12, "-year")
cat("=== IPCW Brier score (apparent) ===\n"); print(round(brier, 3))
write.csv(data.frame(Horizon = rownames(brier), round(brier, 4), check.names = FALSE), file.path(OUT, "table_brier_ipcw.csv"), row.names = FALSE)

# ---- (b) bootstrap calibration slope
slope <- sapply(list(`HAP` = Surv(os_months, event) ~ hap, `HAP+BCLC` = Surv(os_months, event) ~ hap + bclc,
                     `Fixed three-feature model` = Surv(os_months, event) ~ z1 + z2 + z3 + hap + bclc), function(fo) {
  f <- cph(fo, data = d, x = TRUE, y = TRUE); set.seed(42)
  v <- validate(f, method = "boot", B = 500)
  c(Apparent_slope = unname(v["Slope", "index.orig"]), Optimism = unname(v["Slope", "optimism"]), Bootstrap_corrected_slope = unname(v["Slope", "index.corrected"]))
})
cat("\n=== Calibration slope of the linear predictor (rms::validate, B = 500) ===\n"); print(round(t(slope), 2))
write.csv(data.frame(Model = colnames(slope), round(t(slope), 3), row.names = NULL), file.path(OUT, "table_calibration_slope_bootstrap.csv"), row.names = FALSE)
