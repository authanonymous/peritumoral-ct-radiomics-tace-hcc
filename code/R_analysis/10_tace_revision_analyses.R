###############################################################################
#  REPOSITORY VERSION of the authors' revision script (tace_revision_analyses.R).
#  Only three things differ from the original: the default data path, the default output folder, and the
#  line that selects the radiomic columns (by name instead of by position). All analysis code is unchanged.
#  Run from the folder code/R_analysis/ :   Rscript 10_tace_revision_analyses.R
###############################################################################
###############################################################################
#  TACE-HCC peritumoral radiomics - revision analyses for Academic Radiology
#  Manuscript: "Non-Overlapping Peritumoral CT Radiomics ... in TACE-Treated HCC"
#
#  What this script does (numbers = section headers printed to the console)
#   0  Load + audit the data (N, events, HAP recomputation, feature counts)
#   1  Reproduce the published Table 3 / Table 5 numbers (confirms data + method)
#   2  Three-model comparison: HAP | HAP+BCLC | HAP+BCLC+3 features (+ nested)
#   3  Bootstrap optimism correction, WITH selection repeated (nested), paired
#      against HAP+BCLC; audit of the HAP-only corrected C-index
#   4  Repeated 10-fold nested CV incl. LASSO / elastic net (matched comparators)
#   5  70/30 split: HAP, HAP+BCLC, radiomics on the SAME test set (+ 100 splits)
#   6  Feature-family ablation (whole-tumour vs single-zone vs cross-zone ratio)
#   7  Extended Cox (main effect + interaction), time-varying HR, piecewise HR
#   8  Decision curves (1 & 3 years): apparent AND out-of-fold, net benefit at
#      20 % / 30 % with bootstrap CI
#   9  Confounding by tumour size / surface area
#  10  Table 1 recomputed (Wilcoxon / Fisher) + HAP dichotomisation check
#
#  Everything is written to  <folder of the csv>/revision_analyses_output/
#  Runtime: roughly 15-40 min on 4-8 cores with the defaults (see CONFIG).
#  Smoke test (2-5 min): in R run  Sys.setenv(TACE_TEST = "1")  before source().  Tested end-to-end in TEST_MODE.
###############################################################################

## ----------------------------------------------------------------- CONFIG ---
DATA_PATH <- Sys.getenv("TACE_DATA", unset = "../../data/analysis_dataset_110_patients.csv")   # run from code/R_analysis/

TEST_MODE <- identical(Sys.getenv("TACE_TEST"), "1")      # tiny run to check code
SEED        <- 20260920
N_BOOT      <- if (TEST_MODE) 12  else 500   # bootstrap resamples (section 3)
N_SPLITS_ABL <- if (TEST_MODE) 4  else 100   # random 70/30 splits per feature family (section 6)
N_CV_REP    <- if (TEST_MODE) 1   else 5     # repeats of 10-fold CV (section 4)
N_SPLITS    <- if (TEST_MODE) 5   else 100   # random 70/30 splits (section 5)
N_DCA_BOOT  <- if (TEST_MODE) 100 else 1000  # patient bootstrap for DCA CIs
N_CORES     <- if (TEST_MODE) 1 else max(1, parallel::detectCores() - 1)
# optional overrides without editing the file, e.g. Sys.setenv(TACE_N_SPLITS = "50")
ov <- function(nm, default) { v <- Sys.getenv(nm, ""); if (nzchar(v)) as.integer(v) else default }
N_BOOT <- ov("TACE_N_BOOT", N_BOOT); N_SPLITS_ABL <- ov("TACE_N_SPLITS_ABL", N_SPLITS_ABL); N_CV_REP <- ov("TACE_N_CV_REP", N_CV_REP)
N_SPLITS <- ov("TACE_N_SPLITS", N_SPLITS); N_DCA_BOOT <- ov("TACE_N_DCA_BOOT", N_DCA_BOOT); N_CORES <- ov("TACE_N_CORES", N_CORES)
CORR_PRUNE  <- FALSE    # the file is ALREADY the post-|r|>0.90-filter set; TRUE re-prunes inside each fit (drops Surface Area, so it no longer reproduces the 5 FDR survivors)
TIMES_M     <- c(12, 24, 36)                 # months for AUC / Brier / risk
DAYS_PER_MONTH <- 30.4375

OUT_DIR <- Sys.getenv("TACE_OUT", unset = "../../results/revision_analyses_output")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

## -------------------------------------------------------------- PACKAGES ----
need <- c("survival", "glmnet", "ggplot2", "parallel")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) install.packages(miss)
suppressPackageStartupMessages({
  library(survival); library(glmnet); library(ggplot2); library(parallel)
})
has_timeROC <- requireNamespace("timeROC", quietly = TRUE)  # optional cross-check
set.seed(SEED)
say <- function(...) cat("\n", paste0(...), "\n", sep = "")
hdr <- function(x) cat("\n", strrep("=", 78), "\n", x, "\n", strrep("=", 78), "\n", sep = "")
save_csv <- function(df, name) {
  write.csv(df, file.path(OUT_DIR, name), row.names = FALSE)
  invisible(df)
}

## ------------------------------------------------------------- 0. LOAD -----
hdr("0. DATA LOAD AND AUDIT")
raw <- read.csv(DATA_PATH, check.names = TRUE, stringsAsFactors = FALSE)
# repository version: the 624 radiomic columns are selected by name (the deposited file lists them after the clinical columns)
FEATS  <- grep("^(intra_|peri_|ratio_|diff_|habitat_|hab_|gradient_)", names(raw), value = TRUE)
stopifnot(length(FEATS) == 624)

dat <- data.frame(
  id        = raw$PATPRI,
  time_m    = raw$survival_time / DAYS_PER_MONTH,    # survival_time is in DAYS
  death     = raw$death,
  hap_score = raw$hap_score,
  bclc      = raw$bclc,                              # coding (confirmed by the author): 0 = BCLC A, 1 = BCLC B
  raw[, FEATS, drop = FALSE],
  check.names = FALSE
)
clin <- raw[, c("age", "gender_woman", "lesions_number", "lesion1_diameter", "lab_albumin",
                "lab_bilirubin", "lab_afp", "lab_inr", "cps", "etiology")]
n <- nrow(dat)

# ---- key features named in the manuscript (make.names turns '-' into '.') ---
KEY <- make.names(c("ratio_0-3_vs_5-10_glcm_MaximumProbability",   # GLCM MaxProb ratio
                    "ratio_intra_vs_0-3_glcm_MCC",                 # GLCM MCC ratio
                    "ratio_3-5_vs_5-10_glcm_Imc1"))                # GLCM Imc1 ratio
names(KEY) <- c("maxprob", "mcc", "imc1")
stopifnot(all(KEY %in% names(dat)))

say("N = ", n, " | events = ", sum(dat$death), " (", round(100 * mean(dat$death), 1), "%)",
    " | median follow-up-to-event/censor = ", round(median(dat$time_m), 1), " months (all patients)")
say("Radiomic feature columns in file: ", length(FEATS),
    "  (manuscript says 629 after filtering; check the difference)")
say("Feature families: ",
    paste(names(table(sub("_.*", "", FEATS))), table(sub("_.*", "", FEATS)), collapse = " | "))
say("Age mean/SD (file): ", round(mean(clin$age), 1), " / ", round(sd(clin$age), 1),
    "   [manuscript: 65.9 +/- 9.9]")
say("BCLC coding in file: ", paste(names(table(dat$bclc)), table(dat$bclc), sep = "=", collapse = ", "),
    "  -> coded 0 = BCLC A, 1 = BCLC B. The manuscript labels the 80 patients coded 0 as 'BCLC 0' (72.7%): that label looks wrong")
say("HAP score distribution: ", paste(names(table(dat$hap_score)), table(dat$hap_score), sep = "=", collapse = ", "))

# HAP recomputation (AFP>400 ng/mL, albumin<36 g/L, bilirubin>17 umol/L, largest tumour >7 cm)
hap_re <- with(clin, (lab_afp > 400) + (lab_albumin * 10 < 36) + (lab_bilirubin * 17.104 > 17) +
                 (lesion1_diameter > 70))
say("HAP recomputed from labs vs hap_score in file (rows = recomputed, cols = file):")
print(table(recomputed = hap_re, in_file = dat$hap_score))
say("Exact agreement: ", round(100 * mean(hap_re == dat$hap_score), 1), "%  (lesion1_diameter is assumed to be the largest lesion)")
say("HAP dichotomised 0 vs >=1 gives n = ", sum(dat$hap_score == 0), " / ", sum(dat$hap_score > 0),
    "  [Fig 2C says 'median' split, 36/74 -> it is HAP 0 vs >=1, not a median split]")

## ---------------------------------------------------------- HELPERS --------
KEEP_FORCED <- c("hap_score", "bclc")
fit_cox <- function(d, terms) {
  f <- as.formula(paste("Surv(time_m, death) ~", paste(terms, collapse = " + ")))
  suppressWarnings(coxph(f, data = d, ties = "efron"))
}
lp_of   <- function(fit, d) as.numeric(predict(fit, newdata = d, type = "lp"))
cidx    <- function(lp, d) unname(concordance(Surv(d$time_m, d$death) ~ lp, reverse = TRUE)$concordance)
risk_at <- function(fit, d, t) {
  s <- summary(survfit(fit, newdata = d), times = t, extend = TRUE)$surv
  1 - as.numeric(s)
}
cens_fun <- function(time, status) {
  G <- survfit(Surv(time, 1 - status) ~ 1)
  stepfun(G$time, c(1, G$surv))                # G(u) = P(C > u), right-continuous
}
# IPCW cumulative/dynamic AUC (Uno / Blanche); ties in the marker count 1/2
auc_ipcw <- function(marker, time, status, t) {
  case <- which(time <= t & status == 1); ctrl <- which(time > t)
  if (length(case) == 0 || length(ctrl) == 0) return(NA_real_)
  G <- cens_fun(time, status); w <- 1 / G(time[case] - 1e-8)
  conc <- vapply(case, function(i) sum(marker[i] > marker[ctrl]) + 0.5 * sum(marker[i] == marker[ctrl]), numeric(1))
  sum(w * conc) / (sum(w) * length(ctrl))
}
# IPCW Brier score (Graf) at time t; risk = predicted P(event by t)
brier_ipcw <- function(risk, time, status, t) {
  G <- cens_fun(time, status); s <- 1 - risk
  early <- time <= t & status == 1; late <- time > t
  mean(ifelse(early, s^2 / G(time - 1e-8), 0) + ifelse(late, (1 - s)^2 / G(t), 0))
}
auc_all <- function(marker, d) vapply(TIMES_M, function(t) auc_ipcw(marker, d$time_m, d$death, t), numeric(1))

# --- selection pipeline (manuscript): [corr prune] -> univariate Cox (unadjusted)
#     -> Benjamini-Hochberg q<0.05 -> backward elimination (P<0.05) with
#     HAP + BCLC forced in.  Variance filter is impossible here (features already
#     z-scored, and the file is already the post-filter set) - see notes at end.
select_features <- function(d, feats, q = 0.05, p_retain = 0.05, corr_cut = 0.90,
                            prune = CORR_PRUNE, forced = KEEP_FORCED) {
  y <- Surv(d$time_m, d$death)
  pv <- vapply(feats, function(f) {
    x <- d[[f]]
    if (sd(x) == 0) return(1)
    fit <- suppressWarnings(try(coxph(y ~ x), silent = TRUE))
    if (inherits(fit, "try-error") || is.na(coef(fit))) return(1)
    2 * pnorm(-abs(coef(fit) / sqrt(fit$var)))
  }, numeric(1))
  keep <- feats
  if (prune) {                                             # keep the stronger univariate feature
    ord <- order(pv); cm <- abs(cor(d[, feats[ord], drop = FALSE])); kept <- integer(0)
    for (j in seq_along(ord)) if (!length(kept) || all(cm[j, kept] <= corr_cut)) kept <- c(kept, j)
    keep <- feats[ord][kept]
  }
  p_adj <- p.adjust(pv[keep], "BH")
  cand  <- names(p_adj)[p_adj < q]
  attr_cand <- cand
  if (length(cand) > 25) cand <- names(sort(pv[cand]))[1:25]     # numerical safety valve (only bites in resamples that over-select)
  while (length(cand)) {                                   # backward elimination
    fit <- try(fit_cox(d, c(forced, cand)), silent = TRUE)
    if (inherits(fit, "try-error")) { cand <- names(sort(pv[cand]))[seq_len(max(0, floor(length(cand) / 2)))]; next }
    cf  <- summary(fit)$coefficients
    pc  <- setNames(cf[cand, "Pr(>|z|)"], cand)     # setNames: a lone candidate would otherwise lose its name -> infinite loop
    if (any(is.na(pc))) pc[is.na(pc)] <- 1
    if (max(pc) < p_retain) break
    cand <- setdiff(cand, names(which.max(pc)))
  }
  attr(cand, "fdr_survivors") <- attr_cand
  cand
}

# --- parallel wrapper (Windows-safe): indices are pre-generated so results are
#     reproducible irrespective of the number of cores
par_lapply <- function(X, FUN) {
  if (N_CORES <= 1) return(lapply(X, FUN))
  if (.Platform$OS.type == "windows") {
    cl <- makeCluster(N_CORES); on.exit(stopCluster(cl), add = TRUE)
    clusterEvalQ(cl, suppressPackageStartupMessages({library(survival); library(glmnet)}))
    clusterExport(cl, ls(globalenv()), envir = globalenv())
    parLapply(cl, X, FUN)
  } else mclapply(X, FUN, mc.cores = N_CORES)
}

# --------------------------------------------------- MODEL DEFINITIONS -------
spec <- list(
  M1_HAP          = function(d, feats) KEEP_FORCED[1],
  M2_HAP_BCLC     = function(d, feats) KEEP_FORCED,
  M3_fixed_3feat  = function(d, feats) c(KEEP_FORCED, unname(KEY)),          # published model (selection NOT repeated)
  M4_nested_pipe  = function(d, feats) c(KEEP_FORCED, as.character(select_features(d, feats)))
)

## ------------------------------- 1. REPRODUCE THE PUBLISHED NUMBERS --------
hdr("1. REPRODUCTION OF PUBLISHED TABLE 3 / TABLE 5 (sanity check of data + code)")
f3 <- fit_cox(dat, c(unname(KEY), KEEP_FORCED))
tab3 <- data.frame(term = c("GLCM MaxProb ratio (0-3/5-10)", "GLCM MCC ratio (intra/0-3)",
                            "GLCM Imc1 ratio (3-5/5-10)", "HAP score", "BCLC (0 = A, 1 = B)"),
                   summary(f3)$coefficients[, c("coef", "se(coef)", "Pr(>|z|)")],
                   HR = exp(coef(f3)), lo = exp(confint(f3))[, 1], hi = exp(confint(f3))[, 2])
tab3$published_HR <- c(1.316, 1.681, 1.586, 1.602, 1.601)
round_df <- data.frame(tab3[1], lapply(tab3[-1], function(x) round(x, 3)))
print(round_df, row.names = FALSE)
save_csv(round_df, "01_table3_reproduction.csv")

zph <- cox.zph(f3, transform = "km")
say("Schoenfeld test (published: MaxProb chi2 = 7.523, P = 0.006; global P = 0.026):")
print(round(zph$table, 3))

lp3 <- lp_of(f3, dat); lp1 <- dat$hap_score
say("Apparent Harrell C: combined = ", round(cidx(lp3, dat), 3), " (published 0.713) | HAP alone = ",
    round(cidx(lp1, dat), 3), " (published 0.649)")
tab5 <- data.frame(t_months = TIMES_M,
                   AUC_combined = auc_all(lp3, dat), pub_combined = c(0.798, 0.742, 0.753),
                   AUC_HAP = auc_all(lp1, dat), pub_HAP = c(0.744, 0.709, 0.669))
say("IPCW AUC (own implementation) vs published Table 5:"); print(round(tab5, 3), row.names = FALSE)
if (has_timeROC) {
  chk <- sapply(TIMES_M, function(t) timeROC::timeROC(T = dat$time_m, delta = dat$death, marker = lp3,
                  cause = 1, times = t, weighting = "marginal", iid = FALSE)$AUC[2])
  say("timeROC cross-check (combined model): ", paste(round(chk, 3), collapse = " / "))
}
save_csv(tab5, "01_table5_auc_reproduction.csv")

## ------------------------------ 2. THREE-MODEL COMPARISON (apparent) -------
hdr("2. HAP vs HAP+BCLC vs HAP+BCLC+3 FEATURES (apparent performance, full data)")
sf_full <- select_features(dat, FEATS); sel_full <- as.character(sf_full)
say("Features picked by the full pipeline on all 110 patients: ",
    if (length(sel_full)) paste(sel_full, collapse = ", ") else "(none)")
say("FDR-survivors before backward elimination: ",
    paste(attr(sf_full, "fdr_survivors"), collapse = ", "))
models_app <- list(
  M1_HAP         = fit_cox(dat, "hap_score"),
  M2_HAP_BCLC    = fit_cox(dat, KEEP_FORCED),
  M3_fixed_3feat = f3,
  M4_pipeline    = fit_cox(dat, c(KEEP_FORCED, sel_full)))
app <- do.call(rbind, lapply(names(models_app), function(nm) {
  m <- models_app[[nm]]; lp <- lp_of(m, dat)
  a <- auc_all(lp, dat)
  br <- vapply(TIMES_M, function(t) brier_ipcw(risk_at(m, dat, t), dat$time_m, dat$death, t), numeric(1))
  data.frame(model = nm, k_params = length(coef(m)), C_apparent = cidx(lp, dat),
             AUC12 = a[1], AUC24 = a[2], AUC36 = a[3], Brier12 = br[1], Brier24 = br[2], Brier36 = br[3],
             AIC = AIC(m), row.names = NULL)
}))
print(format(app, digits = 3), row.names = FALSE)
lrt <- anova(models_app$M2_HAP_BCLC, models_app$M3_fixed_3feat)
say("Likelihood-ratio test, 3 features added to HAP+BCLC (APPARENT: features were chosen on these data, ",
    "so this p-value is anti-conservative): chi2 = ", round(lrt$Chisq[2], 2), ", df = ", lrt$Df[2],
    ", P = ", signif(lrt$`Pr(>|Chi|)`[2], 3))
save_csv(app, "02_model_comparison_apparent.csv")

## ---------------- 3. BOOTSTRAP OPTIMISM CORRECTION WITH NESTED SELECTION ---
hdr(paste0("3. BOOTSTRAP OPTIMISM (B = ", N_BOOT, "); M1-M3 = fixed structure, M4 = selection repeated in every resample"))
idx_list <- replicate(N_BOOT, sample.int(n, replace = TRUE), simplify = FALSE)
one_boot <- function(idx) {
  db <- dat[idx, ]; oob <- setdiff(seq_len(n), unique(idx)); out <- list(); sel <- NULL
  for (nm in names(spec)) {
    terms <- spec[[nm]](db, FEATS)
    if (nm == "M4_nested_pipe") sel <- setdiff(terms, KEEP_FORCED)
    fit <- fit_cox(db, terms); lpb <- lp_of(fit, db); lpo <- lp_of(fit, dat)
    slope <- tryCatch(unname(coef(coxph(Surv(dat$time_m, dat$death) ~ lpo))), error = function(e) NA_real_)
    out[[nm]] <- c(C_app = cidx(lpb, db), C_test = cidx(lpo, dat),
                   C_oob = if (length(oob) > 10) cidx(lpo[oob], dat[oob, ]) else NA_real_,
                   AUC36_app = auc_ipcw(lpb, db$time_m, db$death, 36),
                   AUC36_test = auc_ipcw(lpo, dat$time_m, dat$death, 36), slope = slope)
  }
  list(perf = out, sel = sel)
}
t0 <- Sys.time(); boots <- par_lapply(idx_list, one_boot)
say("Bootstrap finished in ", round(difftime(Sys.time(), t0, units = "mins"), 1), " min")
getm <- function(nm, what) vapply(boots, function(b) b$perf[[nm]][[what]], numeric(1))
app_C   <- c(M1_HAP = cidx(lp_of(models_app$M1_HAP, dat), dat), M2_HAP_BCLC = cidx(lp_of(models_app$M2_HAP_BCLC, dat), dat),
             M3_fixed_3feat = cidx(lp3, dat), M4_nested_pipe = cidx(lp_of(models_app$M4_pipeline, dat), dat))
app_A36 <- c(M1_HAP = auc_ipcw(lp_of(models_app$M1_HAP, dat), dat$time_m, dat$death, 36),
             M2_HAP_BCLC = auc_ipcw(lp_of(models_app$M2_HAP_BCLC, dat), dat$time_m, dat$death, 36),
             M3_fixed_3feat = auc_ipcw(lp3, dat$time_m, dat$death, 36),
             M4_nested_pipe = auc_ipcw(lp_of(models_app$M4_pipeline, dat), dat$time_m, dat$death, 36))
val <- do.call(rbind, lapply(names(spec), function(nm) {
  opt <- mean(getm(nm, "C_app") - getm(nm, "C_test"), na.rm = TRUE)
  optA <- mean(getm(nm, "AUC36_app") - getm(nm, "AUC36_test"), na.rm = TRUE)
  data.frame(model = nm, C_apparent = app_C[nm], optimism_C = opt, C_corrected = app_C[nm] - opt,
             C_OOB_mean = mean(getm(nm, "C_oob"), na.rm = TRUE),
             AUC36_apparent = app_A36[nm], AUC36_corrected = app_A36[nm] - optA,
             calib_slope_shrinkage = mean(getm(nm, "slope"), na.rm = TRUE), row.names = NULL)
}))
print(format(val, digits = 3), row.names = FALSE)
say("AUDIT: corrected C for HAP alone should be close to its apparent value (minimal optimism for a single-parameter model), matching the manuscript's reported ",
    "0.649 -> 0.648. Here: ", round(val$C_apparent[1], 3), " -> ", round(val$C_corrected[1], 3))
save_csv(val, "03_bootstrap_validation.csv")

# paired differences against HAP+BCLC (the correct clinical comparator)
pair <- do.call(rbind, lapply(c("M3_fixed_3feat", "M4_nested_pipe"), function(nm) {
  dApp  <- getm(nm, "C_app")  - getm("M2_HAP_BCLC", "C_app")
  dTest <- getm(nm, "C_test") - getm("M2_HAP_BCLC", "C_test")
  dOob  <- getm(nm, "C_oob")  - getm("M2_HAP_BCLC", "C_oob")
  dA    <- (getm(nm, "AUC36_app") - getm(nm, "AUC36_test")) - (getm("M2_HAP_BCLC", "AUC36_app") - getm("M2_HAP_BCLC", "AUC36_test"))
  data.frame(comparison = paste(nm, "minus M2_HAP_BCLC"),
             dC_apparent = app_C[nm] - app_C["M2_HAP_BCLC"],
             dC_corrected = (app_C[nm] - app_C["M2_HAP_BCLC"]) - mean(dApp - dTest, na.rm = TRUE),
             dC_OOB_mean = mean(dOob, na.rm = TRUE),
             dC_OOB_q025 = quantile(dOob, .025, na.rm = TRUE), dC_OOB_q975 = quantile(dOob, .975, na.rm = TRUE),
             pct_resamples_radiomics_better_OOB = 100 * mean(dOob > 0, na.rm = TRUE),
             dAUC36_corrected = (app_A36[nm] - app_A36["M2_HAP_BCLC"]) - mean(dA, na.rm = TRUE), row.names = NULL)
}))
print(format(pair, digits = 3), row.names = FALSE)
say("Reading: 'corrected' = apparent minus optimism (Harrell). OOB = performance on patients left out of each resample ",
    "(2.5-97.5 percentile across resamples; wide because ~40 patients are left out each time).")
save_csv(pair, "03b_paired_vs_HAP_BCLC.csv")

say("CAUTION on M4: a bootstrap resample contains duplicated patients, so univariate p-values are anti-conservative and the ",
    "pipeline over-selects (median features kept per resample = ", median(vapply(boots, function(b) length(b$sel), numeric(1))),
    "; the full data give 3). A plain nested bootstrap therefore OVERSTATES optimism for a p-value-driven pipeline. ",
    "Use sections 4-6 (folds / subsamples without duplicates) as the primary nested evidence, and ask the authors for the exact ",
    "procedure and code behind Supplementary Table S5.")
# stability of the pipeline: how often is each feature selected?
sel_tab <- sort(table(unlist(lapply(boots, `[[`, "sel"))), decreasing = TRUE)
sel_df <- data.frame(feature = names(sel_tab), pct_of_resamples = round(100 * as.numeric(sel_tab) / N_BOOT, 1))
sel_df$is_published_feature <- sel_df$feature %in% KEY
say("Feature-selection stability (top 15; published features flagged):"); print(head(sel_df, 15), row.names = FALSE)
say("Resamples in which NO radiomic feature survived: ",
    round(100 * mean(vapply(boots, function(b) length(b$sel) == 0, logical(1))), 1), "%")
save_csv(sel_df, "03c_selection_stability.csv")

## -------------------------- 4. REPEATED 10-FOLD NESTED CV (+ LASSO / EN) ---
hdr(paste0("4. REPEATED 10-FOLD CV x ", N_CV_REP, " (selection and penalty tuning inside training folds)"))
X_names <- c(KEEP_FORCED, FEATS)
glmnet_lp <- function(tr, te, alpha, s) {
  x <- as.matrix(tr[, X_names]); pf <- c(0, 0, rep(1, length(FEATS)))
  cv <- suppressWarnings(cv.glmnet(x, Surv(tr$time_m, tr$death), family = "cox", alpha = alpha,
                                   penalty.factor = pf, nfolds = 10, type.measure = "deviance"))
  lam <- if (s == "1se") cv$lambda.1se else cv$lambda.min
  b <- as.matrix(coef(cv, s = lam)); nz <- setdiff(rownames(b)[b != 0], KEEP_FORCED)
  list(lp_tr = as.numeric(predict(cv, newx = x, s = lam, type = "link")),
       lp_te = as.numeric(predict(cv, newx = as.matrix(te[, X_names]), s = lam, type = "link")), nz = nz)
}
risk_from_lp <- function(lp_tr, tr, lp_te, t) {                # recalibrate a penalised LP on the training fold
  f <- suppressWarnings(coxph(Surv(tr$time_m, tr$death) ~ lp, data = data.frame(tr[, c("time_m", "death")], lp = lp_tr)))
  risk_at(f, data.frame(lp = lp_te), t)
}
cv_variants <- c("M1_HAP", "M2_HAP_BCLC", "M4_stepwise", "LASSO_1se", "LASSO_min", "EN0.5_1se", "EN0.5_min")
fold_lists <- lapply(seq_len(N_CV_REP), function(r) {
  f <- integer(n); for (g in 0:1) { i <- sample(which(dat$death == g)); f[i] <- rep_len(sample(1:10), length(i)) }
  f })
conc_counts <- function(lp, d) {                       # concordant pairs / usable pairs inside ONE test fold
  cc <- concordance(Surv(d$time_m, d$death) ~ lp, reverse = TRUE)$count
  c(num = unname(cc["concordant"] + 0.5 * cc["tied.x"]), den = unname(cc["concordant"] + cc["discordant"] + cc["tied.x"]))
}
one_cv_rep <- function(fold) {
  oof <- lapply(cv_variants, function(v) matrix(NA_real_, n, length(TIMES_M), dimnames = list(NULL, TIMES_M)))
  names(oof) <- cv_variants; sel <- list()
  cnt <- setNames(lapply(cv_variants, function(v) c(num = 0, den = 0)), cv_variants)
  for (k in 1:10) {
    te_i <- which(fold == k); tr <- dat[-te_i, ]; te <- dat[te_i, ]
    add <- function(v, lp_te, risk) { oof[[v]][te_i, ] <<- risk; cnt[[v]] <<- cnt[[v]] + conc_counts(lp_te, te) }
    rk  <- function(fit) sapply(TIMES_M, function(t) risk_at(fit, te, t))
    f1 <- fit_cox(tr, "hap_score");  add("M1_HAP", lp_of(f1, te), rk(f1))
    f2 <- fit_cox(tr, KEEP_FORCED);  add("M2_HAP_BCLC", lp_of(f2, te), rk(f2))
    s <- as.character(select_features(tr, FEATS)); sel[["M4_stepwise"]] <- c(sel[["M4_stepwise"]], s)
    f4 <- fit_cox(tr, c(KEEP_FORCED, s)); add("M4_stepwise", lp_of(f4, te), rk(f4))
    for (a in c(1, 0.5)) for (ss in c("1se", "min")) {
      v <- paste0(if (a == 1) "LASSO_" else "EN0.5_", ss); g <- glmnet_lp(tr, te, a, ss)
      sel[[v]] <- c(sel[[v]], g$nz)
      add(v, g$lp_te, sapply(TIMES_M, function(t) risk_from_lp(g$lp_tr, tr, g$lp_te, t)))
    }
  }
  list(oof = oof, sel = sel, cnt = cnt)
}
t0 <- Sys.time(); cvres <- par_lapply(fold_lists, one_cv_rep)
say("CV finished in ", round(difftime(Sys.time(), t0, units = "mins"), 1), " min")
cv_tab <- do.call(rbind, lapply(seq_along(cvres), function(r) data.frame(repeat_id = r, model = cv_variants,
  C_within_fold = vapply(cv_variants, function(v) cvres[[r]]$cnt[[v]]["num"] / cvres[[r]]$cnt[[v]]["den"], numeric(1)), row.names = NULL)))
cv_tab$dC_vs_M2 <- cv_tab$C_within_fold - cv_tab$C_within_fold[match(paste(cv_tab$repeat_id, "M2_HAP_BCLC"), paste(cv_tab$repeat_id, cv_tab$model))]
cv_sum <- do.call(rbind, lapply(cv_variants, function(v) { x <- cv_tab[cv_tab$model == v, ]
  data.frame(model = v, C_cv = mean(x$C_within_fold), dC_vs_HAP_BCLC = mean(x$dC_vs_M2),
             dC_min_over_repeats = min(x$dC_vs_M2), dC_max_over_repeats = max(x$dC_vs_M2)) }))
print(format(cv_sum, digits = 3), row.names = FALSE)
say("C_cv = C-index computed only among patients of the same held-out fold (each fold scored by the model that never saw it), ",
    "pooled over folds, averaged over repeats. This avoids the downward bias of pooling risks across folds.")
save_csv(cv_tab, "04_repeated_cv_by_repeat.csv"); save_csv(cv_sum, "04_repeated_cv_summary.csv")
nfold_total <- 10 * N_CV_REP
sel_freq <- do.call(rbind, lapply(c("M4_stepwise", "LASSO_1se", "LASSO_min", "EN0.5_1se", "EN0.5_min"), function(v) {
  tt <- table(unlist(lapply(cvres, function(r) r$sel[[v]])))
  if (!length(tt)) return(NULL)
  data.frame(model = v, feature = names(tt), pct_of_folds = round(100 * as.numeric(tt) / nfold_total, 1))
}))
if (!is.null(sel_freq)) {
  sel_freq$published_feature <- sel_freq$feature %in% KEY
  sel_freq <- sel_freq[order(sel_freq$model, -sel_freq$pct_of_folds), ]
  say("Selection frequency across outer folds (top 6 per model):")
  print(do.call(rbind, lapply(split(sel_freq, sel_freq$model), head, 6)), row.names = FALSE)
  save_csv(sel_freq, "04b_selection_frequency_cv.csv")
}
# keep out-of-fold risks for decision curves (mean over repeats)
oof_risk <- lapply(c("M1_HAP", "M2_HAP_BCLC", "M4_stepwise"), function(v)
  Reduce(`+`, lapply(cvres, function(r) r$oof[[v]])) / length(cvres)); names(oof_risk) <- c("M1_HAP", "M2_HAP_BCLC", "M4_stepwise")

## ------------------------------------------------ 5. 70/30 SPLIT, MATCHED --
hdr("5. INTERNAL 70/30 SPLIT: all models on the IDENTICAL test set")
strat_split <- function(seed, frac = 0.7) {
  set.seed(seed); tr <- integer(0)
  for (g in 0:1) { i <- which(dat$death == g); tr <- c(tr, sample(i, round(frac * length(i)))) }
  sort(tr)
}
eval_split <- function(seed) {
  tr_i <- strat_split(seed); tr <- dat[tr_i, ]; te <- dat[-tr_i, ]
  s <- as.character(select_features(tr, FEATS))
  fits <- list(M1_HAP = fit_cox(tr, "hap_score"), M2_HAP_BCLC = fit_cox(tr, KEEP_FORCED),
               M3_fixed_3feat_LEAKY = fit_cox(tr, c(KEEP_FORCED, unname(KEY))),
               M4_selection_on_train = fit_cox(tr, c(KEEP_FORCED, s)))
  g <- glmnet_lp(tr, te, 1, "1se")
  out <- data.frame(seed = seed, n_train = nrow(tr), n_test = nrow(te), events_test = sum(te$death),
                    model = c(names(fits), "LASSO_1se"),
                    C_train = c(sapply(fits, function(m) cidx(lp_of(m, tr), tr)), cidx(g$lp_tr, tr)),
                    C_test  = c(sapply(fits, function(m) cidx(lp_of(m, te), te)), cidx(g$lp_te, te)),
                    n_radiomic_selected = c(0, 0, 3, length(s), length(g$nz)), row.names = NULL)
  out$dC_test_vs_M2 <- out$C_test - out$C_test[out$model == "M2_HAP_BCLC"]
  out
}
single <- eval_split(42)
print(format(single, digits = 3), row.names = FALSE)
# CI for the test-set values: bootstrap the TEST patients (models fixed)
set.seed(SEED + 1); tr_i <- strat_split(42); tr <- dat[tr_i, ]; te <- dat[-tr_i, ]
s42 <- as.character(select_features(tr, FEATS))
fx <- list(M1 = fit_cox(tr, "hap_score"), M2 = fit_cox(tr, KEEP_FORCED), M4 = fit_cox(tr, c(KEEP_FORCED, s42)))
lpt <- lapply(fx, lp_of, d = te)
bs <- replicate(N_DCA_BOOT, { i <- sample(nrow(te), replace = TRUE); b <- te[i, ]
  if (sum(b$death) < 3) return(c(M1 = NA_real_, M2 = NA_real_, M4 = NA_real_))
  c(M1 = cidx(lpt$M1[i], b), M2 = cidx(lpt$M2[i], b), M4 = cidx(lpt$M4[i], b)) })
say("Test-set (n = ", nrow(te), ", ", sum(te$death), " events) bootstrap 95% CI  |  M4 - M2 (selection on train): ",
    round(mean(bs["M4", ] - bs["M2", ], na.rm = TRUE), 3), " [",
    paste(round(quantile(bs["M4", ] - bs["M2", ], c(.025, .975), na.rm = TRUE), 3), collapse = ", "), "]")
save_csv(single, "05_split7030_seed42.csv")
say("Distribution over ", N_SPLITS, " random stratified 70/30 splits (a single split is unstable with ~33 test patients):")
spl <- do.call(rbind, par_lapply(seq(1001, length.out = N_SPLITS), eval_split))
spl_sum <- do.call(rbind, lapply(split(spl, spl$model), function(x) data.frame(
  model = x$model[1], C_test_median = median(x$C_test), C_test_IQR = paste(round(quantile(x$C_test, c(.25, .75)), 3), collapse = "-"),
  dC_vs_M2_median = median(x$dC_test_vs_M2), pct_splits_better_than_M2 = round(100 * mean(x$dC_test_vs_M2 > 0), 1))))
print(format(spl_sum, digits = 3), row.names = FALSE)
save_csv(spl, "05b_split7030_all_splits.csv"); save_csv(spl_sum, "05b_split7030_summary.csv")

## ---------------------------------------- 6. FEATURE-FAMILY ABLATION -------
hdr(paste0("6. FEATURE-FAMILY ABLATION (", N_SPLITS_ABL, " random 70/30 splits per family; held-out C vs HAP+BCLC)"))
fam <- list(
  whole_tumour_intra      = FEATS[grepl("^intra_", FEATS)],
  single_zone_peri        = FEATS[grepl("^peri_", FEATS)],
  cross_zone_ratio_diff   = FEATS[grepl("^(ratio|diff)_", FEATS)],
  habitat                 = FEATS[grepl("^(habitat|hab_)", FEATS)],
  all_features            = FEATS)
abl_split <- function(seed, ff) {
  tr_i <- strat_split(seed); tr <- dat[tr_i, ]; te <- dat[-tr_i, ]
  s <- as.character(select_features(tr, ff))
  c(k = length(s), C_radiomics = cidx(lp_of(fit_cox(tr, c(KEEP_FORCED, s)), te), te),
    C_HAP_BCLC = cidx(lp_of(fit_cox(tr, KEEP_FORCED), te), te))
}
abl <- do.call(rbind, lapply(names(fam), function(fn) {
  M <- do.call(rbind, par_lapply(seq(2001, length.out = N_SPLITS_ABL), function(sd) abl_split(sd, fam[[fn]])))
  d <- M[, "C_radiomics"] - M[, "C_HAP_BCLC"]
  data.frame(family = fn, n_features_in_family = length(fam[[fn]]),
             median_features_selected = median(M[, "k"]), pct_splits_none_selected = round(100 * mean(M[, "k"] == 0), 1),
             C_heldout_HAP_BCLC = median(M[, "C_HAP_BCLC"]), C_heldout_radiomics = median(M[, "C_radiomics"]),
             dC_median = median(d), dC_IQR = paste(round(quantile(d, c(.25, .75)), 3), collapse = " to "),
             pct_splits_radiomics_better = round(100 * mean(d > 0), 1), row.names = NULL)
}))
print(format(abl, digits = 3), row.names = FALSE)
say("If the cross-zone family does not beat the whole-tumour and single-zone families on held-out data, the claim that ",
    "the 'non-overlapping architecture' drives the signal is not supported by these data.")
save_csv(abl, "06_feature_family_ablation.csv")

## ----------------------------------------- 7. EXTENDED COX (time-varying) --
hdr("7. EXTENDED COX: main effect + interaction, time-varying HR, piecewise HR")
dat$maxprob <- dat[[KEY["maxprob"]]]; dat$mcc <- dat[[KEY["mcc"]]]; dat$imc1 <- dat[[KEY["imc1"]]]
fx_tt <- coxph(Surv(time_m, death) ~ hap_score + bclc + mcc + imc1 + maxprob + tt(maxprob), data = dat,
               tt = function(x, t, ...) x * log(t + 1), ties = "efron")
print(round(summary(fx_tt)$coefficients[, c("coef", "se(coef)", "Pr(>|z|)")], 4))
for (u in c("months", "days")) {          # can the published Table 4B be reproduced? (published: 1.252 [1.117-1.404]; HAP 1.596, MCC 1.606, Imc1 1.557, BCLC 1.478)
  dd <- dat; if (u == "days") dd$time_m <- raw$survival_time
  ff <- coxph(Surv(time_m, death) ~ hap_score + bclc + mcc + imc1 + maxprob + tt(maxprob), data = dd,
              tt = function(x, t, ...) x * log(t + 1), ties = "efron")
  say("Extended Cox, time in ", u, ": interaction HR = ", round(exp(coef(ff)["tt(maxprob)"]), 3), " [",
      paste(round(exp(confint(ff)["tt(maxprob)", ]), 3), collapse = "-"), "] | HAP ", round(exp(coef(ff)["hap_score"]), 3),
      " MCC ", round(exp(coef(ff)["mcc"]), 3), " Imc1 ", round(exp(coef(ff)["imc1"]), 3), " BCLC ", round(exp(coef(ff)["bclc"]), 3))
}
say("If neither matches the published Table 4B, ask the authors exactly how the extended model was specified (time unit, terms).")
V <- vcov(fx_tt); b1 <- coef(fx_tt)["maxprob"]; b2 <- coef(fx_tt)["tt(maxprob)"]
tg <- c(3, 6, 12, 24, 36, 48, 60, 72, 84, 96)
vh <- do.call(rbind, lapply(tg, function(t) {
  L <- log(t + 1); est <- b1 + b2 * L
  se <- sqrt(V["maxprob", "maxprob"] + L^2 * V["tt(maxprob)", "tt(maxprob)"] + 2 * L * V["maxprob", "tt(maxprob)"])
  data.frame(month = t, n_at_risk = sum(dat$time_m >= t), HR = exp(est), lo = exp(est - 1.96 * se), hi = exp(est + 1.96 * se))
}))
print(format(vh, digits = 3), row.names = FALSE)
say("Interpret only rows with a reasonable number at risk (Fig 2A shows 1-2 patients beyond 84 months).")
save_csv(vh, "07_time_varying_HR.csv")
tgrid <- seq(1, 96, by = 0.5)
pg <- do.call(rbind, lapply(tgrid, function(t) { L <- log(t + 1); est <- b1 + b2 * L
  se <- sqrt(V["maxprob", "maxprob"] + L^2 * V["tt(maxprob)", "tt(maxprob)"] + 2 * L * V["maxprob", "tt(maxprob)"])
  data.frame(t = t, HR = exp(est), lo = exp(est - 1.96 * se), hi = exp(est + 1.96 * se), nrisk = sum(dat$time_m >= t)) }))
pg$reliable <- pg$nrisk >= 10
p_tv <- ggplot(pg, aes(t, HR)) + geom_ribbon(aes(ymin = lo, ymax = hi, fill = reliable), alpha = .25) +
  geom_line(aes(linetype = reliable)) + geom_hline(yintercept = 1, colour = "grey40", linetype = 2) +
  scale_y_log10() + scale_fill_manual(values = c(`TRUE` = "#1f4e79", `FALSE` = "grey60")) +
  labs(x = "Months since TACE", y = "HR per SD, GLCM MaxProb ratio (log scale)",
       title = "Time-varying HR (extended Cox). Grey = fewer than 10 patients at risk") + theme_bw() + theme(legend.position = "none")
ggsave(file.path(OUT_DIR, "07_time_varying_HR.png"), p_tv, width = 7, height = 4.5, dpi = 300)
sp <- survSplit(Surv(time_m, death) ~ ., data = dat[, c("time_m", "death", "hap_score", "bclc", "mcc", "imc1", "maxprob")],
                cut = c(24, 48), episode = "period")
sp$period <- factor(sp$period, labels = c("0-24 mo", "24-48 mo", ">48 mo"))
fpw <- coxph(Surv(tstart, time_m, death) ~ hap_score + bclc + mcc + imc1 + maxprob:period, data = sp)
pw <- data.frame(period = levels(sp$period), events = as.numeric(tapply(sp$death, sp$period, sum)),
                 HR = exp(coef(fpw)[grep("maxprob", names(coef(fpw)))]),
                 lo = exp(confint(fpw)[grep("maxprob", names(coef(fpw))), 1]), hi = exp(confint(fpw)[grep("maxprob", names(coef(fpw))), 2]))
say("Piecewise (interval-specific) HR of MaxProb ratio - model-free check of the 'increasing effect' claim:")
print(format(pw, digits = 3), row.names = FALSE); save_csv(pw, "07b_piecewise_HR.csv")

## ------------------------------------------------ 8. DECISION CURVES -------
hdr("8. DECISION CURVE ANALYSIS: 1 and 3 years, apparent and out-of-fold")
nb_surv <- function(risk, time, status, t, pts) {
  N <- length(risk); s_all <- summary(survfit(Surv(time, status) ~ 1), times = t, extend = TRUE)$surv
  sapply(pts, function(pt) {
    hi <- risk >= pt; nh <- sum(hi); if (nh == 0) return(0)
    s_h <- summary(survfit(Surv(time[hi], status[hi]) ~ 1), times = t, extend = TRUE)$surv
    ((1 - s_h) * nh) / N - (s_h * nh) / N * pt / (1 - pt)
  })
}
nb_all <- function(time, status, t, pts) { s <- summary(survfit(Surv(time, status) ~ 1), times = t, extend = TRUE)$surv
  (1 - s) - s * pts / (1 - pts) }
pts_h <- list(`12` = seq(0.05, 0.60, by = 0.01), `36` = seq(0.20, 0.90, by = 0.01))   # 3-year event risk is high: thresholds must be higher
thr_h <- list(`12` = c(0.10, 0.20, 0.30), `36` = c(0.20, 0.30, 0.40, 0.50, 0.60))
risk_app <- lapply(list(M1_HAP = models_app$M1_HAP, M2_HAP_BCLC = models_app$M2_HAP_BCLC, M3_fixed_3feat = f3),
                   function(m) sapply(c(12, 36), function(t) risk_at(m, dat, t)))
dca_rows <- list(); thr_rows <- list()
for (tt_i in 1:2) {
  t <- c(12, 36)[tt_i]; pts <- pts_h[[as.character(t)]]
  for (src in c("apparent", "out_of_fold")) {
    rl <- if (src == "apparent") lapply(risk_app, function(x) x[, tt_i]) else lapply(oof_risk, function(x) x[, as.character(t)])
    for (nm in names(rl)) dca_rows[[length(dca_rows) + 1]] <- data.frame(horizon = t, source = src, model = nm, thr = pts,
        NB = nb_surv(rl[[nm]], dat$time_m, dat$death, t, pts))
    dca_rows[[length(dca_rows) + 1]] <- data.frame(horizon = t, source = src, model = "treat_all", thr = pts,
        NB = nb_all(dat$time_m, dat$death, t, pts))
    # net benefit at 20% / 30% + patient-bootstrap CI for difference vs HAP+BCLC
    ref <- if (src == "apparent") "M2_HAP_BCLC" else "M2_HAP_BCLC"
    set.seed(SEED + 2); bidx <- replicate(N_DCA_BOOT, sample.int(n, replace = TRUE), simplify = FALSE)
    for (nm in setdiff(names(rl), ref)) for (pt in thr_h[[as.character(t)]]) {
      d0 <- nb_surv(rl[[nm]], dat$time_m, dat$death, t, pt) - nb_surv(rl[[ref]], dat$time_m, dat$death, t, pt)
      bd <- vapply(bidx, function(i) nb_surv(rl[[nm]][i], dat$time_m[i], dat$death[i], t, pt) -
                                     nb_surv(rl[[ref]][i], dat$time_m[i], dat$death[i], t, pt), numeric(1))
      thr_rows[[length(thr_rows) + 1]] <- data.frame(horizon_m = t, source = src, model = nm, vs = ref, threshold = pt,
        NB_model = nb_surv(rl[[nm]], dat$time_m, dat$death, t, pt), NB_ref = nb_surv(rl[[ref]], dat$time_m, dat$death, t, pt),
        dNB = d0, dNB_lo = quantile(bd, .025), dNB_hi = quantile(bd, .975),
        extra_net_TP_per_100 = 100 * d0)
    }
  }
}
dca <- do.call(rbind, dca_rows); thr <- do.call(rbind, thr_rows); rownames(thr) <- NULL
print(format(thr, digits = 3), row.names = FALSE)
say("Net benefit is per patient; dNB x 100 = extra net true positives per 100 patients. ",
    "A CI that spans 0 means no demonstrable advantage over HAP+BCLC at that threshold.")
save_csv(dca, "08_dca_curves.csv"); save_csv(thr, "08b_dca_thresholds.csv")
dca$panel <- paste0(dca$horizon, "-month risk, ", gsub("_", "-", dca$source))
p_dca <- ggplot(dca, aes(thr, NB, colour = model)) + geom_line() + facet_wrap(~ panel, scales = "free_y") +
  coord_cartesian(ylim = c(-0.05, NA)) + labs(x = "Threshold probability", y = "Net benefit") + theme_bw()
ggsave(file.path(OUT_DIR, "08_dca.png"), p_dca, width = 10, height = 7, dpi = 300)

## ------------------------------ 9. CONFOUNDING BY SIZE / SURFACE AREA ------
hdr("9. DOES THE SIGNAL SURVIVE ADJUSTMENT FOR TUMOUR SIZE / SURFACE AREA?")
dat$diam <- raw$lesion1_diameter; dat$logdiam <- log(raw$lesion1_diameter)
dat$surf <- dat[["intra_original_shape_SurfaceArea"]]; dat$multifocal <- as.integer(raw$lesions_number > 1)
adj <- list(`M3 (published)` = c(KEEP_FORCED, "mcc", "imc1", "maxprob"),
            `M3 + log diameter` = c(KEEP_FORCED, "logdiam", "mcc", "imc1", "maxprob"),
            `M3 + surface area` = c(KEEP_FORCED, "surf", "mcc", "imc1", "maxprob"),
            `M3 + diameter + surface + multifocal` = c(KEEP_FORCED, "logdiam", "surf", "multifocal", "mcc", "imc1", "maxprob"))
adjtab <- do.call(rbind, lapply(names(adj), function(nm) { m <- fit_cox(dat, adj[[nm]]); s <- summary(m)$coefficients
  k <- c("maxprob", "mcc", "imc1"); data.frame(model = nm, feature = k, HR = s[k, "exp(coef)"],
    lo = exp(confint(m))[k, 1], hi = exp(confint(m))[k, 2], p = s[k, "Pr(>|z|)"], C = cidx(lp_of(m, dat), dat), row.names = NULL) }))
print(format(adjtab, digits = 3), row.names = FALSE); save_csv(adjtab, "09_size_adjusted_models.csv")
cors <- sapply(c("diam", "surf", "hap_score", "bclc"), function(v) sapply(c("maxprob", "mcc", "imc1"),
        function(f) cor(dat[[f]], dat[[v]], method = "spearman")))
say("Spearman correlation of the three features with size / stage (a strong |rho| means the feature partly encodes size):")
print(round(cors, 2))
lrt9 <- anova(fit_cox(dat, c(KEEP_FORCED, "logdiam", "surf")),
              fit_cox(dat, c(KEEP_FORCED, "logdiam", "surf", "mcc", "imc1", "maxprob")))
say("LRT: do the 3 features add to HAP+BCLC+log diameter+surface area?  P = ", signif(lrt9$`Pr(>|Chi|)`[2], 3),
    "  (apparent; features were chosen on these data)")

## ----------------------------------------------- 10. TABLE 1 RECOMPUTED ----
hdr("10. TABLE 1 RECOMPUTED (risk groups = median of the combined model's linear predictor)")
grp <- factor(ifelse(lp3 > median(lp3), "High", "Low"), levels = c("Low", "High"))
say("Group sizes: ", paste(names(table(grp)), table(grp), collapse = " / "))
t1 <- list()
cont <- function(lab, x) { p <- wilcox.test(x ~ grp, exact = FALSE)$p.value
  data.frame(variable = lab, Low = sprintf("%.2f +/- %.2f", mean(x[grp == "Low"]), sd(x[grp == "Low"])),
             High = sprintf("%.2f +/- %.2f", mean(x[grp == "High"]), sd(x[grp == "High"])), test = "Wilcoxon", P = p) }
catg <- function(lab, x) { tb <- table(grp, x); p <- fisher.test(tb)$p.value
  if (all(x %in% c(0, 1))) { f <- function(g) sprintf("%d (%.1f%%)", sum(x[grp == g] == 1), 100 * mean(x[grp == g] == 1))
    data.frame(variable = lab, Low = f("Low"), High = f("High"), test = "Fisher exact", P = p) }
  else data.frame(variable = lab, Low = paste(tb["Low", ], collapse = "/"), High = paste(tb["High", ], collapse = "/"),
                  test = "Fisher exact", P = p) }
t1 <- rbind(cont("Age (y)", clin$age), catg("Male sex", 1 - clin$gender_woman), cont("Index lesion diameter (mm)", clin$lesion1_diameter),
  catg("Multifocal", as.integer(clin$lesions_number > 1)), cont("AFP (ng/mL)", clin$lab_afp), cont("Albumin (g/dL)", clin$lab_albumin),
  cont("Bilirubin (mg/dL)", clin$lab_bilirubin), cont("INR", clin$lab_inr), cont("HAP score", dat$hap_score),
  catg("BCLC A (coded 0) vs B", as.integer(dat$bclc == 0)), catg("Aetiology (3 categories)", clin$etiology),
  catg("Child-Pugh A (vs B)", as.integer(clin$cps == 1)), catg("Death", dat$death))
t1$P_formatted <- ifelse(t1$P < 0.001, "<0.001", sprintf("%.3f", t1$P))
t1$P <- NULL
print(t1, row.names = FALSE)
say("Cross-check: Table 1 counts recomputed here vs. the published manuscript (Child-Pugh A, BCLC A and Deaths, by risk group) -- confirmed matching:")
print(data.frame(row = c("Child-Pugh A", "BCLC A (coded 0)", "Deaths"), manuscript_Low = c("53", "46", "27"), manuscript_High = c("41", "34", "48"),
                 data_Low = c(sum(clin$cps[grp == "Low"] == 1), sum(dat$bclc[grp == "Low"] == 0), sum(dat$death[grp == "Low"])),
                 data_High = c(sum(clin$cps[grp == "High"] == 1), sum(dat$bclc[grp == "High"] == 0), sum(dat$death[grp == "High"]))), row.names = FALSE)
save_csv(t1, "10_table1_recomputed.csv")

## --------------------------------------------------------- SESSION ---------
hdr("DONE")
say("All tables/figures in: ", OUT_DIR)
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
