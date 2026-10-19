# =============================================================================
# 00_common.R - shared helper functions (sourced by scripts 01-04). Do not run alone.
# All scripts: run from this folder (code/R_analysis/); input ../../data/analysis_dataset_110_patients.csv
# Packages: survival (all), glmnet (03), rms is not needed here.
# =============================================================================
suppressPackageStartupMessages(library(survival))

DATA <- "../../data/analysis_dataset_110_patients.csv"
OUT  <- "../../results"
if (!file.exists(DATA)) stop("Run the scripts from code/R_analysis/ - cannot find: ", DATA)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# Fixed three-feature model (chosen once on the full cohort; see 01_univariate_fdr_stepwise.R)
FIXED3 <- c("ratio_intra_vs_0-3_glcm_MCC",
            "ratio_0-3_vs_5-10_glcm_MaximumProbability",
            "ratio_3-5_vs_5-10_glcm_Imc1")
SHORT <- c("ratio_intra_vs_0-3_glcm_MCC" = "MCC ratio (intra/0-3)",
           "ratio_0-3_vs_5-10_glcm_MaximumProbability" = "Maximum Probability ratio (0-3/5-10)",
           "ratio_3-5_vs_5-10_glcm_Imc1" = "Imc1 ratio (3-5/5-10)")

load_data <- function() {
  d <- read.csv(DATA, check.names = FALSE, stringsAsFactors = FALSE)
  d$os_months <- d$survival_time / 30.44          # days -> months
  d$event     <- as.integer(d$death)
  d$hap       <- as.numeric(d$hap_score)
  d$bclc      <- as.numeric(d$bclc)                # 0 = A, 1 = B
  d
}
feature_names <- function(d) {
  f <- grep("^(intra_|peri_|ratio_|diff_|habitat_|hab_|gradient_)", names(d), value = TRUE)
  stopifnot(length(f) == 624)                      # 624 features after variance + correlation filtering
  f
}
bt <- function(x) sprintf("`%s`", x)               # back-quote non-syntactic names for formulas

# Univariate Cox for every feature (features are already z-scored: HR per SD), Wald P, Benjamini-Hochberg q.
# NOTE: all quantities are computed from the data frame passed in (never from objects outside it),
# so the function can safely be used on bootstrap samples, folds and splits.
uni_screen <- function(d, feats) {
  tt <- d$os_months; ee <- d$event
  res <- t(vapply(feats, function(f) {
    x <- d[[f]]
    fit <- tryCatch(suppressWarnings(coxph(Surv(tt, ee) ~ x)), error = function(e) NULL)
    if (is.null(fit) || !is.finite(coef(fit)) || !is.finite(fit$var[1, 1])) return(c(NA_real_, NA_real_, 1))
    s <- summary(fit)$coefficients
    c(s[1, "coef"], s[1, "se(coef)"], s[1, "Pr(>|z|)"])
  }, numeric(3)))
  out <- data.frame(Feature = feats, coef = res[, 1], se = res[, 2], p_value = res[, 3], stringsAsFactors = FALSE)
  out$HR <- exp(out$coef); out$CI_low <- exp(out$coef - 1.96 * out$se); out$CI_high <- exp(out$coef + 1.96 * out$se)
  out$FDR_q <- p.adjust(out$p_value, method = "BH")
  out$Significant_FDR <- ifelse(out$FDR_q < 0.05, "Yes", "No")
  rownames(out) <- NULL
  out
}

# Selection pipeline: univariate Cox + BH (q < q_cut), then backward elimination of the radiomic candidates
# (P < p_cut) with HAP score and BCLC stage always in the model. If nothing passes, the model is HAP + BCLC.
select_pipeline <- function(d, feats, q_cut = 0.05, p_cut = 0.05) {
  scr <- uni_screen(d, feats)
  sel <- scr$Feature[!is.na(scr$FDR_q) & scr$FDR_q < q_cut]
  repeat {
    terms <- c(if (length(sel)) bt(sel), "hap", "bclc")
    fit <- tryCatch(suppressWarnings(coxph(reformulate(terms, response = "Surv(os_months, event)"), data = d)),
                    error = function(e) NULL)
    if (is.null(fit)) {                               # degenerate fit: fall back to the clinical model
      sel <- character(0)
      fit <- coxph(Surv(os_months, event) ~ hap + bclc, data = d); break
    }
    if (!length(sel)) break
    pv <- summary(fit)$coefficients[, "Pr(>|z|)"]
    names(pv) <- gsub("`", "", names(pv), fixed = TRUE)
    pv <- pv[sel]; pv[is.na(pv)] <- 1
    if (all(pv < p_cut)) break
    sel <- setdiff(sel, names(which.max(pv)))
  }
  list(features = sel, fit = fit, screen = scr)
}

fit_fixed3 <- function(d) coxph(reformulate(c(bt(FIXED3), "hap", "bclc"), response = "Surv(os_months, event)"), data = d)

score_lp <- function(fit, newdata) as.numeric(predict(fit, newdata = newdata, type = "lp"))

# Harrell C: larger linear predictor = higher risk = shorter survival
cindex <- function(lp, d) {
  as.numeric(concordance(Surv(d$os_months, d$event) ~ lp, reverse = TRUE)$concordance)
}

# Stratified (by event status) K-fold assignment
strat_folds <- function(event, K, seed) {
  set.seed(seed); id <- integer(length(event))
  for (s in unique(event)) { ix <- sample(which(event == s)); id[ix] <- rep_len(sample(seq_len(K)), length(ix)) }
  id
}
# Stratified train/test split (train fraction = frac)
strat_split <- function(event, frac, seed) {
  set.seed(seed); train <- logical(length(event))
  for (s in unique(event)) { ix <- which(event == s); train[sample(ix, round(frac * length(ix)))] <- TRUE }
  train
}
