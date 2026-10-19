# =============================================================================
# 02_models_ph_extended_cox.R
# (a) adjustment for tumour burden (log diameter, surface area, multifocality)         -> Supplementary Table S4
# (b) Schoenfeld proportional-hazards tests                                            -> PH section of Results
# (c) extended Cox model: Maximum Probability ratio x log(t + 1), t in months         -> Table 3B, Fig 4
# (d) time-specific HRs (3, 12, 24, 48 months) and interval-specific HRs               -> Supplementary Table S5
# Expected (manuscript): PH violation for MaxProb (chi2 7.523, P = 0.006; global P = 0.026); extended model main effect
#   HR 0.61 (0.31-1.22), interaction HR 1.64 (1.13-2.39, P = 0.009); time-specific HR 1.22 / 2.18 / 3.02 / 4.22;
#   interval HRs 1.29 (0-24 mo), 3.44 (24-48 mo), 2.24 (>48 mo).
# =============================================================================
source("00_common.R")
d <- load_data(); MP <- "ratio_0-3_vs_5-10_glcm_MaximumProbability"
d$z1 <- d[[MP]]; d$z2 <- d[["ratio_intra_vs_0-3_glcm_MCC"]]; d$z3 <- d[["ratio_3-5_vs_5-10_glcm_Imc1"]]
d$log_diam <- log(d$lesion1_diameter); d$surface <- d[["intra_original_shape_SurfaceArea"]]; d$multifocal <- as.integer(d$lesions_number > 1)

base <- Surv(os_months, event) ~ z1 + z2 + z3 + hap + bclc
cat("=== (a) Adjustment for tumour burden: HR (95% CI) for the three features ===\n")
adj <- list("Fixed three-feature model" = base,
            "+ log tumour diameter" = update(base, . ~ . + log_diam),
            "+ surface area" = update(base, . ~ . + surface),
            "+ log diameter + surface area + multifocality" = update(base, . ~ . + log_diam + surface + multifocal))
rows <- do.call(rbind, lapply(names(adj), function(nm) {
  f <- coxph(adj[[nm]], data = d); s <- summary(f)$coefficients; ci <- exp(confint(f))
  v <- c("z1", "z2", "z3")
  data.frame(Model = nm, Feature = c("Maximum Probability ratio", "MCC ratio", "Imc1 ratio"), HR = s[v, "exp(coef)"],
             CI_low = ci[v, 1], CI_high = ci[v, 2], P = s[v, "Pr(>|z|)"], row.names = NULL)
}))
print(rows, digits = 3, row.names = FALSE); write.csv(rows, file.path(OUT, "table_S4_adjusted_models.csv"), row.names = FALSE)

cat("\n=== (b) Scaled Schoenfeld residual tests (fixed three-feature model) ===\n")
f0 <- coxph(base, data = d); zph <- cox.zph(f0); print(zph)
write.csv(data.frame(Term = rownames(zph$table), zph$table, row.names = NULL), file.path(OUT, "table_PH_schoenfeld.csv"), row.names = FALSE)

cat("\n=== (c) Extended Cox model with z1 x log(t + 1) (t in months) ===\n")
ext <- coxph(Surv(os_months, event) ~ z1 + z2 + z3 + hap + bclc + tt(z1), data = d,
             tt = function(x, t, ...) x * log(t + 1))
s <- summary(ext)$coefficients; ci <- exp(confint(ext)); print(cbind(HR = s[, "exp(coef)"], ci, P = s[, "Pr(>|z|)"]), digits = 3)

cat("\n=== (d) Time-specific HR per SD of the Maximum Probability ratio ===\n")
V <- vcov(ext); b1 <- coef(ext)["z1"]; b2 <- coef(ext)["tt(z1)"]
tsp <- do.call(rbind, lapply(c(3, 12, 24, 48), function(m) {
  l <- log(m + 1); est <- b1 + b2 * l; se <- sqrt(V["z1", "z1"] + l^2 * V["tt(z1)", "tt(z1)"] + 2 * l * V["z1", "tt(z1)"])
  data.frame(Month = m, HR = exp(est), CI_low = exp(est - 1.96 * se), CI_high = exp(est + 1.96 * se))
}))
print(tsp, digits = 3, row.names = FALSE)

cat("\nInterval-specific HR (0-24, 24-48, >48 months; other covariates common):\n")
ds <- survSplit(Surv(os_months, event) ~ ., data = d, cut = c(24, 48), episode = "tgroup")
ds$tgroup <- factor(ds$tgroup, labels = c("0-24", "24-48", ">48"))
fi <- coxph(Surv(tstart, os_months, event) ~ z2 + z3 + hap + bclc + z1:strata(tgroup), data = ds)
si <- summary(fi)$coefficients; cii <- exp(confint(fi)); ix <- grep("^z1", rownames(si))
ev <- tapply(ds$event, ds$tgroup, sum)
itv <- data.frame(Interval = levels(ds$tgroup), HR = si[ix, "exp(coef)"], CI_low = cii[ix, 1], CI_high = cii[ix, 2], Events = as.integer(ev), row.names = NULL)
print(itv, digits = 3, row.names = FALSE)
write.csv(rbind(cbind(Type = "time-specific", tsp), cbind(Type = "interval-specific", setNames(itv[, 1:4], names(tsp)))),
          file.path(OUT, "table_S5_time_varying.csv"), row.names = FALSE)
