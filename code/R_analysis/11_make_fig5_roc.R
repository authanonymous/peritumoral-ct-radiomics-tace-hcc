# =============================================================================
# 11_make_fig5_roc.R - Fig 5: time-dependent ROC curves (inverse-probability-of-censoring weighted, Kaplan-Meier censoring weights,
# i.e. the "marginal" estimator of timeROC) at 1, 2 and 3 years for the fixed three-feature model and HAP alone.
# Prints the AUCs (expected 0.798 / 0.742 / 0.753 and 0.744 / 0.709 / 0.669 = Table 4). Run from code/R_analysis/.
# =============================================================================
suppressPackageStartupMessages({library(survival); library(ggplot2)})
DATA <- "../../data/analysis_dataset_110_patients.csv"
if (!file.exists(DATA)) stop("Run from code/R_analysis/ - cannot find: ", DATA)
d <- read.csv(DATA, check.names = FALSE)
d$tm <- d$survival_time / 30.44; d$ev <- d$death
d$z1 <- d[["ratio_0-3_vs_5-10_glcm_MaximumProbability"]]; d$z2 <- d[["ratio_intra_vs_0-3_glcm_MCC"]]; d$z3 <- d[["ratio_3-5_vs_5-10_glcm_Imc1"]]
lp_fix <- predict(coxph(Surv(tm, ev) ~ z1 + z2 + z3 + hap_score + bclc, data = d), type = "lp")
lp_hap <- predict(coxph(Surv(tm, ev) ~ hap_score, data = d), type = "lp")
# KM estimate of the censoring survival function, left limit G(t-)
G_left <- function(x) { cs <- sort(unique(d$tm[d$ev == 0])); g <- 1
  for (c in cs[cs < x]) g <- g * (1 - sum(d$tm == c & d$ev == 0) / sum(d$tm >= c)); g }
roc_ipcw <- function(marker, h) {
  case <- d$tm <= h & d$ev == 1; ctrl <- d$tm > h
  w <- vapply(which(case), function(i) 1 / G_left(d$tm[i]), 0); w <- w / sum(w)
  th <- sort(unique(marker), decreasing = TRUE)
  tpr <- c(0, vapply(th, function(c) sum(w * (marker[case] >= c)), 0)); fpr <- c(0, vapply(th, function(c) mean(marker[ctrl] >= c), 0))
  # AUC with ties counted 1/2 (trapezoid over threshold steps)
  list(curve = data.frame(fpr = fpr, tpr = tpr), auc = sum(diff(fpr) * (head(tpr, -1) + tail(tpr, -1)) / 2))
}
rows <- list(); ann <- list()
for (y in 1:3) { a <- roc_ipcw(lp_fix, 12 * y); b <- roc_ipcw(lp_hap, 12 * y)
  cat(sprintf("%d-year AUC: fixed three-feature %.3f | HAP %.3f\n", y, a$auc, b$auc))
  lab <- paste0(y, "-Year ROC")
  rows[[length(rows) + 1]] <- cbind(a$curve, Model = "Fixed three-feature model", Panel = lab)
  rows[[length(rows) + 1]] <- cbind(b$curve, Model = "HAP alone", Panel = lab)
  ann[[y]] <- data.frame(Panel = lab, label = sprintf("Fixed three-feature model AUC=%.3f\nHAP AUC=%.3f", a$auc, b$auc)) }
df <- do.call(rbind, rows); an <- do.call(rbind, ann)
p <- ggplot(df, aes(fpr, tpr, colour = Model, linetype = Model)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey60", linetype = "dotted") +
  geom_path(linewidth = 1.2) +
  geom_text(data = an, aes(x = 0.99, y = 0.02, label = label), inherit.aes = FALSE, hjust = 1, vjust = 0, size = 2.6, colour = "grey25", lineheight = 1.05) +
  facet_wrap(~Panel, nrow = 1) +
  scale_colour_manual(values = c("Fixed three-feature model" = "#1B4F72", "HAP alone" = "#E74C3C")) +
  scale_linetype_manual(values = c("Fixed three-feature model" = "solid", "HAP alone" = "dashed")) +
  labs(x = "1 \u2013 Specificity", y = "Sensitivity", title = "Time-Dependent ROC Curves \u2014 Fixed Three-Feature Model vs HAP") +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = TRUE) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, size = 12), strip.background = element_blank(), strip.text = element_text(face = "bold", hjust = 0, size = 11),
        legend.position = "bottom", legend.key.width = unit(1.6, "cm"), panel.spacing = unit(1.2, "lines"))
dir.create("../../results", showWarnings = FALSE)
ggsave("../../results/Fig5_roc_fixed_three_feature_model.png", p, width = 1432 / 150, height = 1015 / 150, dpi = 150, bg = "white")
