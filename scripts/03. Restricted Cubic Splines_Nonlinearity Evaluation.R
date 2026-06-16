# =============================================================================
# PROJECT: Clinical Prediction Model — Occult Distant Metastasis in Gastric Cancer
# STEP 3:  Restricted Cubic Splines — Nonlinearity Evaluation
# AUTHOR:  Nhan Tran
# DATE:    160626
# DEPENDENCY: Step 1 and Step 2 must be run first.
#             Required objects: df, model_glm, model_lrm, dd,
#                               cindex_apparent, cindex_corrected,
#                               slope_corrected, brier_apparent
# =============================================================================


# =============================================================================
# SECTION 0: Environment & Prerequisites
# =============================================================================

set.seed(2025)

library(rms)
library(ggplot2)
library(dplyr)
library(gt)
library(pROC)

# Bootstrap iterations
B_val <- 2000

# Safety checks
if (!exists("model_lrm")) stop("model_lrm not found. Run Steps 1-2 first.")
if (!exists("df"))        stop("df not found. Run Steps 1-2 first.")

# Confirm datadist
options(datadist = "dd")

cat("=== STEP 3 PREREQUISITES CHECK ===\n")
cat("n =", nrow(df), "| Events =", sum(df$Metastasis_TNM), "\n")
cat("Model A corrected C-index (Step 2 benchmark): ~0.719\n")
cat("Model A corrected calibration slope (Step 2): ~0.831\n\n")


# =============================================================================
# SECTION 1 — PART 1: Exploratory Functional Form Plots
#
# Purpose: Visual inspection of NLR and PLR vs outcome
# Method:  Smoothed scatter (loess) on binary outcome
#          Binned observed probability plot
# These are EXPLORATORY only — do not use alone to conclude nonlinearity
# =============================================================================

cat("=== PART 1: EXPLORATORY FUNCTIONAL FORM PLOTS ===\n")

# --- Helper: binned probability plot ---
# Divides predictor into quantile bins, plots observed event rate per bin
plot_binned_outcome <- function(predictor_vec,
                                outcome_vec,
                                predictor_name,
                                n_bins = 8) {
  
  bin_df <- data.frame(
    x = predictor_vec,
    y = as.numeric(outcome_vec)
  ) %>%
    mutate(
      bin = cut(x,
                breaks   = quantile(x, probs = seq(0, 1, length.out = n_bins + 1)),
                include.lowest = TRUE,
                labels   = FALSE)
    ) %>%
    group_by(bin) %>%
    summarise(
      mean_x    = mean(x),
      obs_prob  = mean(y),
      n         = n(),
      se        = sqrt(obs_prob * (1 - obs_prob) / n),
      .groups   = "drop"
    )
  
  ggplot(bin_df, aes(x = mean_x, y = obs_prob)) +
    geom_point(size = 3, color = "#2C3E8C") +
    geom_errorbar(aes(ymin = obs_prob - 1.96 * se,
                      ymax = obs_prob + 1.96 * se),
                  width = 0, color = "#2C3E8C", alpha = 0.6) +
    geom_smooth(method = "loess", se = TRUE,
                color  = "#E07B39", fill = "#E07B39",
                alpha  = 0.15, linewidth = 1) +
    scale_y_continuous(limits = c(0, 1),
                       labels = scales::percent_format(accuracy = 1)) +
    labs(
      title    = paste0("Observed Metastasis Rate by ", predictor_name),
      subtitle = paste0("Points = binned observed rate (n_bins = ", n_bins,
                        "); Orange = loess smooth"),
      x        = predictor_name,
      y        = "Observed Probability of Metastasis"
    ) +
    theme_classic(base_size = 13) +
    theme(plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(color = "gray40"))
}

# Generate plots
p_nlr_explore <- plot_binned_outcome(df$NLR, df$Metastasis_TNM, "NLR")
p_plr_explore <- plot_binned_outcome(df$PLR, df$Metastasis_TNM, "PLR")

print(p_nlr_explore)
print(p_plr_explore)

ggsave("figure4a_NLR_functional_form.png", p_nlr_explore,
       width = 6, height = 4.5, dpi = 300)
ggsave("figure4b_PLR_functional_form.png", p_plr_explore,
       width = 6, height = 4.5, dpi = 300)

cat("Exploratory figures saved.\n")
cat("Visual inspection: note any U-shape, plateau, or threshold patterns.\n\n")


# =============================================================================
# SECTION 2 — PART 2: Fit RCS Model (Model B)
#
# Model B: Metastasis_TNM ~ Tumor_mlocation3 + rcs(NLR, 4) + rcs(PLR, 4)
# 4 knots placed at default rms quantile positions (5th, 35th, 65th, 95th)
# x = TRUE, y = TRUE required for validate() and calibrate() in Part 6
# =============================================================================

cat("=== PART 2: FITTING RCS MODEL (MODEL B) ===\n")

model_rcs <- lrm(
  Metastasis_TNM ~ Tumor_mlocation3 + rcs(NLR, 4) + rcs(PLR, 4),
  data = df,
  x    = TRUE,
  y    = TRUE
)

cat("\n--- MODEL B SUMMARY ---\n")
print(model_rcs)

# Extract apparent C-index from Model B (lrm gives C directly)
cindex_B_apparent <- model_rcs$stats["C"]

cat(sprintf("\nModel B apparent C-index: %.4f\n", cindex_B_apparent))

# AIC for both models
# Model A: directly from glm object
aic_A_check <- AIC(model_glm)

# Model B: lrm stores -2LL in stats; AIC = -2LL + 2*(df + 1)
# +1 accounts for the intercept term
deviance_A <- -2 * as.numeric(logLik(model_glm))
deviance_B <- model_rcs$stats["-2 LL"]
df_B       <- model_rcs$stats["d.f."]
aic_B      <- deviance_B + 2 * (df_B + 1)

cat(sprintf("AIC Model A (linear):  %.2f\n", aic_A_check))
cat(sprintf("AIC Model B (RCS):     %.2f\n", aic_B))
cat(sprintf("Delta AIC (B - A):     %.2f\n", aic_B - aic_A_check))
cat("(Negative delta AIC favors Model B; >2 units = meaningful difference)\n")

# NOTE: rms::lrm does not export AIC directly via AIC()
# Use logLik-based manual calculation:
# AIC = -2 * logLik + 2 * df
loglik_A <- as.numeric(logLik(model_glm))
df_A     <- length(coef(model_glm))
aic_A    <- -2 * loglik_A + 2 * df_A

loglik_B <- model_rcs$stats["Model L.R."] / (-2) +
  model_rcs$stats["-2 LL"]/ (-2)
# Safer: extract directly from lrm stats
# lrm stores: stats["-2 LL"] = deviance of fitted model
deviance_B <- model_rcs$stats["-2 LL"]
df_B       <- model_rcs$stats["d.f."]     # model df used
aic_B      <- deviance_B + 2 * (df_B + 1) # +1 for intercept

deviance_A <- -2 * loglik_A
aic_A_check <- deviance_A + 2 * df_A

cat(sprintf("\nAIC Model A (linear):  %.2f\n", aic_A_check))
cat(sprintf("AIC Model B (RCS):     %.2f\n", aic_B))
cat(sprintf("Delta AIC (B - A):     %.2f\n", aic_B - aic_A_check))
cat("(Negative delta AIC favors Model B; >2 units = meaningful difference)\n")


# =============================================================================
# SECTION 3 — PART 3: ANOVA Nonlinearity Tests
#
# anova(lrm_model) partitions the LR chi-square into:
#   - Overall association (all terms for that variable)
#   - Nonlinear component (the added spline df beyond linear)
# p(nonlinear) < 0.05 → evidence of nonlinearity
# =============================================================================

cat("\n=== PART 3: ANOVA NONLINEARITY TESTS ===\n")

anova_rcs <- anova(model_rcs)
print(anova_rcs)

# Extract NLR nonlinearity p-value
# anova() for rms lrm returns rows for each term + "TOTAL"
# Row names follow pattern: "rcs(NLR, 4) Nonlinear"
anova_df <- as.data.frame(anova_rcs)
cat("\n--- Formatted ANOVA Table ---\n")
print(anova_df)

# Attempt to extract nonlinear p-values by row name pattern
nlr_rows <- grep("NLR", rownames(anova_df), value = TRUE)
plr_rows <- grep("PLR", rownames(anova_df), value = TRUE)

cat("\nNLR-related ANOVA rows:\n"); print(anova_df[nlr_rows, ])
cat("\nPLR-related ANOVA rows:\n"); print(anova_df[plr_rows, ])

# Manual interpretation guide
cat("\n--- NONLINEARITY INTERPRETATION ---\n")
cat("For each predictor, locate the 'Nonlinear' row:\n")
cat("  p < 0.05 → Evidence of nonlinear relationship\n")
cat("  p ≥ 0.05 → Linear term likely sufficient\n")
cat("  Note: p-value is necessary but NOT sufficient for adopting RCS model\n\n")


# =============================================================================
# SECTION 4 — PART 4: Spline Effect Plots
#
# Uses rms::Predict() to generate predicted log-odds across NLR/PLR range
# Other predictors held at median/reference
# Converted to OR scale for interpretability
# =============================================================================

cat("=== PART 4: SPLINE EFFECT PLOTS ===\n")

# --- NLR spline plot ---
pred_nlr <- Predict(
  model_rcs,
  NLR  = seq(min(df$NLR), max(df$NLR), length.out = 200),
  fun  = exp,          # Convert log-odds to OR
  conf.int = 0.95
)

pred_nlr_df <- as.data.frame(pred_nlr)

p_nlr_spline <- ggplot(pred_nlr_df, aes(x = NLR, y = yhat)) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              fill = "#2C3E8C", alpha = 0.15) +
  geom_line(color = "#2C3E8C", linewidth = 1.2) +
  geom_hline(yintercept = 1, linetype = "dashed",
             color = "gray50") +
  # Rug plot showing data distribution
  geom_rug(data = df, aes(x = NLR, y = NULL),
           sides = "b", alpha = 0.4, color = "gray40") +
  scale_y_log10() +    # Log scale for OR is conventional
  labs(
    title    = "Effect of NLR on Metastasis Risk (RCS, 4 knots)",
    subtitle = "Other predictors held at reference/median",
    x        = "NLR (Neutrophil-to-Lymphocyte Ratio)",
    y        = "Odds Ratio (log scale)",
    caption  = "Shaded area = 95% CI; Rug = observed data distribution"
  ) +
  theme_classic(base_size = 13) +
  theme(plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(color = "gray40"))

print(p_nlr_spline)
ggsave("figure5a_NLR_spline_effect.png", p_nlr_spline,
       width = 6.5, height = 5, dpi = 300)

# --- PLR spline plot ---
pred_plr <- Predict(
  model_rcs,
  PLR  = seq(min(df$PLR), max(df$PLR), length.out = 200),
  fun  = exp,
  conf.int = 0.95
)

pred_plr_df <- as.data.frame(pred_plr)

p_plr_spline <- ggplot(pred_plr_df, aes(x = PLR, y = yhat)) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              fill  = "#B5543B", alpha = 0.15) +
  geom_line(color = "#B5543B", linewidth = 1.2) +
  geom_hline(yintercept = 1, linetype = "dashed",
             color = "gray50") +
  geom_rug(data = df, aes(x = PLR, y = NULL),
           sides = "b", alpha = 0.4, color = "gray40") +
  scale_y_log10() +
  labs(
    title    = "Effect of PLR on Metastasis Risk (RCS, 4 knots)",
    subtitle = "Other predictors held at reference/median",
    x        = "PLR (Platelet-to-Lymphocyte Ratio)",
    y        = "Odds Ratio (log scale)",
    caption  = "Shaded area = 95% CI; Rug = observed data distribution"
  ) +
  theme_classic(base_size = 13) +
  theme(plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(color = "gray40"))

print(p_plr_spline)
ggsave("figure5b_PLR_spline_effect.png", p_plr_spline,
       width = 6.5, height = 5, dpi = 300)

cat("Spline effect figures saved.\n")


# =============================================================================
# SECTION 5 — PART 5: Model A vs Model B Apparent Comparison
#
# Compares on: C-index, AUC, AIC, LR test, apparent Brier
# Bootstrap-corrected comparison in Part 6
# =============================================================================

cat("\n=== PART 5: MODEL A vs MODEL B APPARENT COMPARISON ===\n")

# Model A apparent metrics (from Step 1)
pred_probs_A   <- predict(model_glm, type = "response")
roc_A          <- roc(df$Metastasis_TNM, pred_probs_A,
                      quiet = TRUE, direction = "<")
auc_A_apparent <- round(as.numeric(auc(roc_A)), 4)
brier_A        <- round(mean((pred_probs_A - df$Metastasis_TNM)^2), 4)
cindex_A_app   <- round(as.numeric(auc_A_apparent), 4)  # AUC = C-index

# Model B apparent metrics
pred_probs_B   <- predict(model_rcs, type = "fitted")
roc_B          <- roc(df$Metastasis_TNM, pred_probs_B,
                      quiet = TRUE, direction = "<")
auc_B_apparent <- round(as.numeric(auc(roc_B)), 4)
brier_B        <- round(mean((pred_probs_B - df$Metastasis_TNM)^2), 4)
cindex_B_app   <- round(as.numeric(auc_B_apparent), 4)

# Likelihood ratio test: Model A nested in Model B?
# Model A is nested in Model B (RCS reduces to linear when nonlinear df = 0)
# LR statistic = deviance_A - deviance_B; df = df_B - df_A
# LR test will be computed in Part 7 after deviance_B is confirmed

# Build comparison table
model_comparison <- data.frame(
  Metric      = c("C-index (Apparent)",
                  "AUC (Apparent)",
                  "AIC",
                  "Brier Score (Apparent)",
                  "LR test p-value"),
  Model_A_Linear = c(
    cindex_A_app,
    auc_A_apparent,
    round(aic_A_check, 2),
    brier_A,
    "—"
  ),
  Model_B_RCS    = c(
    cindex_B_app,
    auc_B_apparent,
    round(aic_B, 2),
    brier_B,
    round(lr_p, 4)
  )
)

print(model_comparison)

# gt table
comparison_gt2 <- model_comparison %>%
  gt() %>%
  tab_header(
    title    = "Table 3. Model A (Linear) vs Model B (RCS) — Apparent Performance",
    subtitle = "Full dataset, n = 114"
  ) %>%
  cols_label(
    Metric         = "Metric",
    Model_A_Linear = "Model A (Linear)",
    Model_B_RCS    = "Model B (RCS, 4 knots)"
  ) %>%
  tab_footnote(
    footnote = paste0(
      "AUC = C-index for binary logistic regression. ",
      "AIC: lower = better. ",
      "Brier Score: lower = better. ",
      "LR test: likelihood ratio test of Model B vs Model A (nested). ",
      "RCS = Restricted Cubic Splines with 4 knots (NLR and PLR)."
    )
  )

print(comparison_gt2)


# =============================================================================
# SECTION 6 — PART 6: Bootstrap Validation of Model B
#
# Identical procedure to Step 2, applied to model_rcs
# Allows direct comparison of optimism-corrected metrics between models
# =============================================================================

cat("\n=== PART 6: BOOTSTRAP VALIDATION OF MODEL B ===\n")
cat("Running validate() on Model B with B =", B_val, "iterations...\n")

set.seed(2025)

val_rcs <- validate(
  model_rcs,
  method = "boot",
  B      = B_val
)

print(val_rcs)

# Extract metrics
dxy_B_apparent  <- val_rcs["Dxy", "index.orig"]
dxy_B_corrected <- val_rcs["Dxy", "index.corrected"]
dxy_B_optimism  <- val_rcs["Dxy", "optimism"]

cindex_B_corr   <- dxy_B_corrected / 2 + 0.5
cindex_B_optim  <- dxy_B_optimism  / 2

slope_B_apparent  <- val_rcs["Slope", "index.orig"]
slope_B_corrected <- val_rcs["Slope", "index.corrected"]
slope_B_optimism  <- val_rcs["Slope", "optimism"]

cat(sprintf("\nModel B | C-index: Apparent=%.3f | Optimism=%.3f | Corrected=%.3f\n",
            cindex_B_app, cindex_B_optim, cindex_B_corr))
cat(sprintf("Model B | Cal Slope: Apparent=%.3f | Optimism=%.3f | Corrected=%.3f\n",
            slope_B_apparent, slope_B_optimism, slope_B_corrected))

# --- Head-to-head bootstrap comparison table ---
# Step 2 values — update these with your actual Step 2 results
cindex_A_corrected_step2  <- 0.719   # REPLACE with your actual Step 2 value
slope_A_corrected_step2   <- 0.831   # REPLACE with your actual Step 2 value

head_to_head <- data.frame(
  Metric = c(
    "Apparent C-index",
    "Optimism (C-index)",
    "Corrected C-index",
    "Apparent Cal. Slope",
    "Optimism (Slope)",
    "Corrected Cal. Slope"
  ),
  Model_A = c(
    round(cindex_A_app,               3),
    round(cindex_A_app -
            cindex_A_corrected_step2, 3),   # back-calculated optimism
    round(cindex_A_corrected_step2,   3),
    round(slope_apparent,             3),
    round(slope_apparent -
            slope_A_corrected_step2,  3),
    round(slope_A_corrected_step2,    3)
  ),
  Model_B = c(
    round(cindex_B_app,      3),
    round(cindex_B_optim,    3),
    round(cindex_B_corr,     3),
    round(slope_B_apparent,  3),
    round(slope_B_optimism,  3),
    round(slope_B_corrected, 3)
  )
)

head_to_head_gt <- head_to_head %>%
  gt() %>%
  tab_header(
    title    = "Table 4. Bootstrap Validation: Model A vs Model B",
    subtitle = paste0("B = ", B_val, " bootstrap iterations")
  ) %>%
  cols_label(
    Metric  = "Metric",
    Model_A = "Model A (Linear)",
    Model_B = "Model B (RCS)"
  ) %>%
  tab_footnote(
    footnote = paste0(
      "Model A values from Step 2. ",
      "Corrected = Apparent − Optimism. ",
      "Higher corrected C-index = better discrimination. ",
      "Corrected calibration slope closer to 1.0 = better calibration."
    )
  )

print(head_to_head)
print(head_to_head_gt)


# =============================================================================
# SECTION 7 — PART 7: Model Selection Recommendation
#
# Decision criteria:
#   Adopt RCS (Model B) ONLY if ALL of:
#     1. Nonlinear ANOVA p < 0.05 for at least one predictor
#     2. Corrected C-index gain > 0.01-0.02 over Model A
#     3. Corrected calibration slope not substantially worse
#     4. AIC favors Model B (lower)
#   Otherwise: retain Model A
# =============================================================================

cat("\n=== PART 7: MODEL SELECTION RECOMMENDATION ===\n\n")

# --- Likelihood Ratio Test ---
# Model A deviance from glm; Model B deviance from lrm stats
deviance_A <- -2 * as.numeric(logLik(model_glm))
deviance_B <- as.numeric(model_rcs$stats["-2 LL"])
df_A       <- length(coef(model_glm))
df_B       <- as.numeric(model_rcs$stats["d.f."]) + 1  # +1 for intercept

lr_stat <- deviance_A - deviance_B
lr_df   <- df_B - df_A
lr_p    <- ifelse(lr_df > 0,
                  pchisq(lr_stat, df = lr_df, lower.tail = FALSE),
                  NA_real_)

cat(sprintf("LR test: Chi2 = %.3f, df = %d, p = %.4f\n", lr_stat, lr_df, lr_p))

# --- Recommendation logic with NA guards ---
delta_cindex <- cindex_B_corr - cindex_A_corrected_step2
delta_slope  <- slope_B_corrected - slope_A_corrected_step2
delta_aic    <- aic_B - aic_A_check

cat(sprintf("Delta Corrected C-index (B - A): %+.3f\n", delta_cindex))
cat(sprintf("Delta Corrected Cal. Slope (B - A): %+.3f\n", delta_slope))
cat(sprintf("Delta AIC (B - A): %+.2f\n", delta_aic))

# Safe NA handling before logical tests
lr_p_safe  <- ifelse(is.na(lr_p),  1, lr_p)
aic_safe   <- ifelse(is.na(delta_aic), 0, delta_aic)

adopt_rcs <- (lr_p_safe < 0.05) & (delta_cindex > 0.01) & (aic_safe < 0)

if (adopt_rcs) {
  cat("RECOMMENDATION: ADOPT MODEL B (RCS)\n\n")
  cat("  - ANOVA nonlinearity test significant (p < 0.05)\n")
  cat("  - Corrected C-index improved by >0.01\n")
  cat("  - AIC favors Model B\n")
} else {
  cat("RECOMMENDATION: RETAIN MODEL A (LINEAR)\n\n")
  if (!is.na(lr_p) && lr_p >= 0.05) {
    cat("  - No significant nonlinearity detected (LR test p >= 0.05)\n")
  }
  if (delta_cindex <= 0.01) {
    cat(sprintf("  - Corrected C-index gain negligible (delta = %+.3f)\n", delta_cindex))
  }
  if (!is.na(delta_aic) && delta_aic >= 0) {
    cat(sprintf("  - AIC does not favor Model B (delta AIC = %+.2f)\n", delta_aic))
  }
  cat("  - Added complexity of RCS not justified by data\n")
  cat("  - Model A proceeds as the primary model for Steps 4-10\n")
}

cat("\n  NOTE: Final decision rests with the Strategist AI and Product Owner.\n")
cat("  This recommendation is based on pre-specified decision criteria.\n")


# =============================================================================
# SECTION 8: Step 3 Summary
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("STEP 3 COMPLETE — RCS EVALUATION SUMMARY\n")
cat("================================================================\n")
cat(sprintf("Model A (Linear) — Corrected C-index: %.3f | Cal Slope: %.3f\n",
            cindex_A_corrected_step2, slope_A_corrected_step2))
cat(sprintf("Model B (RCS)    — Corrected C-index: %.3f | Cal Slope: %.3f\n",
            cindex_B_corr, slope_B_corrected))
cat(sprintf("Delta C-index (corrected): %+.3f\n", delta_cindex))
cat(sprintf("Delta AIC: %+.2f\n", delta_aic))
cat(sprintf("LR test p: %.4f\n", lr_p))
cat("---\n")
cat("Figures:  figure4a_NLR_functional_form.png\n")
cat("          figure4b_PLR_functional_form.png\n")
cat("          figure5a_NLR_spline_effect.png\n")
cat("          figure5b_PLR_spline_effect.png\n")
cat("Tables:   Table 3 (apparent comparison)\n")
cat("          Table 4 (bootstrap head-to-head)\n")
cat("================================================================\n")
cat("STOPPING — Awaiting review before Step 4 (Decision Curve Analysis).\n")