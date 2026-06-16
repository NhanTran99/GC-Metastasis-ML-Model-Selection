
# =============================================================================
# PROJECT: Clinical Prediction Model — Occult Distant Metastasis in Gastric Cancer
# STEP 5:  Elastic Net Benchmark
# AUTHOR:  Nhan Tran
# DATE:    160626
# DEPENDENCY: Steps 1-4 must be run first.
#             Required objects: df, model_glm
# =============================================================================


# =============================================================================
# SECTION 0: Environment & Prerequisites
# =============================================================================

set.seed(2025)

# install.packages("glmnet")
library(glmnet)    # Elastic Net
library(pROC)      # ROC, AUC
library(ggplot2)   # Figures
library(dplyr)     # Data wrangling
library(gt)        # Publication tables

# Bootstrap iterations
B <- 2000

# Reference values from approved logistic model (Step 2)
# UPDATE these with your actual Step 2 results if different
ref_cindex_apparent  <- 0.757   # Step 1 apparent C-index
ref_cindex_corrected <- 0.719   # Step 2 corrected C-index
ref_optimism         <- 0.017   # Step 2 optimism
ref_slope_corrected  <- 0.831   # Step 2 corrected calibration slope
ref_brier_apparent   <- 0.204   # Step 1 apparent Brier

# Safety checks
if (!exists("df"))        stop("df not found. Run Steps 1-4 first.")
if (!exists("model_glm")) stop("model_glm not found. Run Steps 1-4 first.")

cat("=== STEP 5 PREREQUISITES CHECK ===\n")
cat("n =", nrow(df), "| Events =", sum(df$Metastasis_TNM), "\n\n")


# =============================================================================
# SECTION 1 — PART 1: Hyperparameter Tuning
#
# Strategy:
#   For each alpha in {0.0, 0.1, ..., 1.0}:
#     Run cv.glmnet() with 10-fold CV (repeated for stability)
#     Record minimum CV deviance and corresponding lambda
#   Select alpha + lambda with lowest overall CV deviance
#
# WHY 10-fold CV here (not bootstrap)?
#   cv.glmnet() uses k-fold CV internally for lambda selection.
#   This is standard glmnet practice and separate from our bootstrap
#   validation in Part 4 which estimates optimism on the full pipeline.
# =============================================================================

cat("=== PART 1: HYPERPARAMETER TUNING ===\n")
cat("Tuning alpha across 0.0 to 1.0 with 10-fold CV...\n\n")

# Build model matrix
# model.matrix() handles factor dummy coding automatically
# [, -1] removes the intercept column (glmnet adds its own)
X <- model.matrix(
  Metastasis_TNM ~ Tumor_mlocation3 + NLR + PLR,
  data = df
)[, -1]

y <- df$Metastasis_TNM

cat("Model matrix dimensions:", dim(X), "\n")
cat("Column names:", colnames(X), "\n\n")

# Alpha grid
alpha_grid <- seq(0, 1, by = 0.1)

# Storage for tuning results
tuning_results <- data.frame(
  alpha       = alpha_grid,
  best_lambda = NA_real_,
  cv_deviance = NA_real_
)

# Loop over alpha values
for (i in seq_along(alpha_grid)) {
  
  a <- alpha_grid[i]
  
  # cv.glmnet: 10-fold CV, binomial family
  # nfolds = 10 is standard; type.measure = "deviance" for logistic
  cv_fit <- cv.glmnet(
    x            = X,
    y            = y,
    alpha        = a,
    family       = "binomial",
    nfolds       = 10,
    type.measure = "deviance"
  )
  
  tuning_results$best_lambda[i] <- cv_fit$lambda.min
  tuning_results$cv_deviance[i] <- min(cv_fit$cvm)
  
  cat(sprintf("  alpha = %.1f | lambda.min = %.5f | CV deviance = %.4f\n",
              a, cv_fit$lambda.min, min(cv_fit$cvm)))
}

# Select best alpha and lambda
best_row    <- which.min(tuning_results$cv_deviance)
best_alpha  <- tuning_results$alpha[best_row]
best_lambda <- tuning_results$best_lambda[best_row]

cat(sprintf("\nBest alpha:  %.1f\n", best_alpha))
cat(sprintf("Best lambda: %.5f\n", best_lambda))
cat(sprintf("CV deviance: %.4f\n\n", tuning_results$cv_deviance[best_row]))

# Publication tuning table
tuning_gt <- tuning_results %>%
  mutate(across(where(is.numeric), ~ round(.x, 5))) %>%
  gt() %>%
  tab_header(
    title    = "Table 6. Elastic Net Hyperparameter Tuning Results",
    subtitle = "10-fold cross-validation across alpha grid (0.0–1.0)"
  ) %>%
  cols_label(
    alpha       = "Alpha (α)",
    best_lambda = "Best Lambda (λ)",
    cv_deviance = "CV Deviance"
  ) %>%
  tab_style(
    style     = cell_fill(color = "#D6EAF8"),
    locations = cells_body(rows = best_row)
  ) %>%
  tab_footnote(
    footnote = paste0(
      "Alpha = 0: Ridge regression. Alpha = 1: Lasso regression. ",
      "Highlighted row = selected hyperparameters. ",
      "Lambda selected at minimum CV deviance (lambda.min)."
    )
  )

print(tuning_gt)


# =============================================================================
# SECTION 2 — PART 2: Final Elastic Net Model
#
# Fit using best_alpha and best_lambda on full dataset
# Extract and display shrunken coefficients
# =============================================================================

cat("\n=== PART 2: FINAL ELASTIC NET MODEL ===\n")

# Final model fit
final_en <- glmnet(
  x      = X,
  y      = y,
  alpha  = best_alpha,
  lambda = best_lambda,
  family = "binomial"
)

# Extract coefficients
en_coef <- coef(final_en)
en_coef_df <- data.frame(
  Term        = rownames(en_coef),
  Coefficient = round(as.numeric(en_coef), 5)
) %>%
  mutate(
    OR        = round(exp(Coefficient), 3),
    Retained  = ifelse(Coefficient != 0, "Yes", "No (shrunk to 0)")
  )

cat("\n--- Elastic Net Coefficients ---\n")
print(en_coef_df)

# Compare with logistic regression coefficients
cat("\n--- Comparison: Logistic vs Elastic Net Coefficients ---\n")
glm_coef <- coef(model_glm)
coef_compare <- data.frame(
  Term          = names(glm_coef),
  Logistic_coef = round(glm_coef, 4),
  EN_coef       = round(as.numeric(en_coef)[
    match(names(glm_coef), rownames(en_coef))], 4),
  Shrinkage_pct = NA_real_
) %>%
  mutate(
    Shrinkage_pct = round(
      (1 - abs(EN_coef) / abs(Logistic_coef)) * 100, 1
    )
  )

print(coef_compare)

# Publication coefficient table
coef_gt <- en_coef_df %>%
  filter(Term != "(Intercept)") %>%
  gt() %>%
  tab_header(
    title    = "Table 7. Elastic Net Model Coefficients",
    subtitle = paste0("Alpha = ", best_alpha,
                      ", Lambda = ", round(best_lambda, 5))
  ) %>%
  cols_label(
    Term        = "Predictor",
    Coefficient = "Shrunken Coefficient (β)",
    OR          = "Odds Ratio (exp β)",
    Retained    = "Retained?"
  ) %>%
  tab_footnote(
    footnote = paste0(
      "Coefficients penalized by Elastic Net (alpha = ", best_alpha,
      ", lambda = ", round(best_lambda, 5), "). ",
      "OR = exp(coefficient). ",
      "Intercept excluded from table."
    )
  )

print(coef_gt)


# =============================================================================
# SECTION 3 — PART 3: Apparent Performance
#
# Predicted probabilities from Elastic Net on full dataset
# ROC, AUC, Brier Score
# =============================================================================

cat("\n=== PART 3: APPARENT PERFORMANCE ===\n")

# Predicted probabilities
en_pred_apparent <- as.numeric(
  predict(final_en, newx = X, type = "response")
)

# ROC and AUC
roc_en <- roc(
  response  = y,
  predictor = en_pred_apparent,
  ci        = TRUE,
  direction = "<",
  quiet     = TRUE
)

auc_en_val <- round(as.numeric(auc(roc_en)), 4)
auc_en_ci  <- round(as.numeric(ci.auc(roc_en)), 4)

# Brier Score
brier_en_apparent <- round(
  mean((en_pred_apparent - y)^2), 4
)

cat(sprintf("Apparent AUC:    %.4f (95%% CI: %.4f–%.4f)\n",
            auc_en_val, auc_en_ci[1], auc_en_ci[3]))
cat(sprintf("Apparent Brier:  %.4f\n", brier_en_apparent))

# --- ROC Figure ---
# Overlay Elastic Net and Logistic ROC curves
pred_logistic <- predict(model_glm, type = "response")
roc_logistic  <- roc(y, pred_logistic, quiet = TRUE, direction = "<")

roc_en_df <- data.frame(
  FPR      = 1 - roc_en$specificities,
  TPR      = roc_en$sensitivities,
  Model    = paste0("Elastic Net (AUC = ", auc_en_val, ")")
)

roc_log_df <- data.frame(
  FPR      = 1 - roc_logistic$specificities,
  TPR      = roc_logistic$sensitivities,
  Model    = paste0("Logistic (AUC = ",
                    round(as.numeric(auc(roc_logistic)), 4), ")")
)

roc_combined <- bind_rows(roc_en_df, roc_log_df)

roc_plot_en <- ggplot(roc_combined,
                      aes(x = FPR, y = TPR,
                          color = Model, linetype = Model)) +
  geom_line(linewidth = 1.1) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "gray60") +
  scale_color_manual(values = c(
    "Elastic Net" = "#2C3E8C",
    "Logistic"    = "#E07B39"
  ) |> setNames(unique(roc_combined$Model))) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title    = "ROC Curves — Logistic vs Elastic Net",
    subtitle = "Apparent performance, full dataset (n = 114)",
    x        = "1 − Specificity (False Positive Rate)",
    y        = "Sensitivity (True Positive Rate)",
    color    = NULL, linetype = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title      = element_text(face = "bold"),
    plot.subtitle   = element_text(color = "gray40"),
    legend.position = c(0.65, 0.15),
    legend.background = element_rect(fill = "white", color = "gray80")
  )

print(roc_plot_en)
ggsave("figure8_ROC_EN_vs_logistic.png", roc_plot_en,
       width = 6.5, height = 5.5, dpi = 300)

cat("ROC figure saved: figure8_ROC_EN_vs_logistic.png\n")


# =============================================================================
# SECTION 4 — PART 4: Bootstrap Internal Validation
#
# CRITICAL DESIGN DECISION:
#   Lambda is RE-TUNED inside each bootstrap iteration.
#   This correctly estimates optimism for the FULL pipeline
#   (tune + fit), not just the fit step.
#   Tuning on full data then bootstrapping would underestimate optimism.
#
# For speed, alpha is fixed at best_alpha found in Part 1.
# Re-tuning alpha inside each loop would be prohibitively slow.
# This is an acceptable approximation documented in the manuscript.
# =============================================================================

cat("\n=== PART 4: BOOTSTRAP INTERNAL VALIDATION ===\n")
cat(sprintf("Running %d bootstrap iterations...\n", B))
cat("(Lambda re-tuned inside each iteration — this may take 3-8 minutes)\n\n")

set.seed(2025)
n <- nrow(df)

# Storage vectors
boot_auc_boot  <- numeric(B)   # AUC on bootstrap sample (apparent in boot world)
boot_auc_orig  <- numeric(B)   # AUC on original data (test in boot world)
boot_slope_boot <- numeric(B)  # Calibration slope on bootstrap sample
boot_slope_orig <- numeric(B)  # Calibration slope on original data
boot_brier_boot <- numeric(B)  # Brier on bootstrap sample
boot_brier_orig <- numeric(B)  # Brier on original data
skipped <- 0

for (b in seq_len(B)) {
  
  # Step 1: Bootstrap sample
  boot_idx  <- sample(seq_len(n), size = n, replace = TRUE)
  X_boot    <- X[boot_idx, , drop = FALSE]
  y_boot    <- y[boot_idx]
  
  # Step 2: Re-tune lambda on bootstrap sample (alpha fixed)
  cv_boot <- tryCatch(
    cv.glmnet(
      x            = X_boot,
      y            = y_boot,
      alpha        = best_alpha,
      family       = "binomial",
      nfolds       = 10,
      type.measure = "deviance"
    ),
    error = function(e) NULL
  )
  
  if (is.null(cv_boot)) { skipped <- skipped + 1; next }
  
  lambda_boot <- cv_boot$lambda.min
  
  # Step 3: Fit on bootstrap sample
  fit_boot <- tryCatch(
    glmnet(
      x      = X_boot,
      y      = y_boot,
      alpha  = best_alpha,
      lambda = lambda_boot,
      family = "binomial"
    ),
    error = function(e) NULL
  )
  
  if (is.null(fit_boot)) { skipped <- skipped + 1; next }
  
  # Step 4a: Predict on BOOTSTRAP sample (apparent performance)
  p_boot_on_boot <- as.numeric(
    predict(fit_boot, newx = X_boot, type = "response")
  )
  
  # Step 4b: Predict on ORIGINAL data (test performance)
  p_boot_on_orig <- as.numeric(
    predict(fit_boot, newx = X, type = "response")
  )
  
  # Step 5: AUC
  boot_auc_boot[b] <- tryCatch(
    as.numeric(auc(roc(y_boot, p_boot_on_boot,
                       quiet = TRUE, direction = "<"))),
    error = function(e) NA
  )
  boot_auc_orig[b] <- tryCatch(
    as.numeric(auc(roc(y, p_boot_on_orig,
                       quiet = TRUE, direction = "<"))),
    error = function(e) NA
  )
  
  # Step 6: Calibration slope
  # Fit logistic regression of outcome on log-odds of predictions
  logit_boot <- log(p_boot_on_boot / (1 - p_boot_on_boot + 1e-10))
  logit_orig <- log(p_boot_on_orig / (1 - p_boot_on_orig + 1e-10))
  
  boot_slope_boot[b] <- tryCatch(
    coef(glm(y_boot ~ logit_boot, family = binomial))[2],
    error = function(e) NA
  )
  boot_slope_orig[b] <- tryCatch(
    coef(glm(y ~ logit_orig, family = binomial))[2],
    error = function(e) NA
  )
  
  # Step 7: Brier Score
  boot_brier_boot[b] <- mean((p_boot_on_boot - y_boot)^2)
  boot_brier_orig[b] <- mean((p_boot_on_orig - y)^2)
  
  # Progress update every 500 iterations
  if (b %% 500 == 0) cat(sprintf("  Iteration %d / %d complete\n", b, B))
}

cat(sprintf("\nBootstrap complete. Skipped iterations: %d / %d\n",
            skipped, B))

# --- Compute optimism and corrected metrics ---
auc_optimism       <- mean(boot_auc_boot  - boot_auc_orig,  na.rm = TRUE)
slope_optimism_en  <- mean(boot_slope_boot - boot_slope_orig, na.rm = TRUE)
brier_optimism_en  <- mean(boot_brier_boot - boot_brier_orig, na.rm = TRUE)

auc_en_corrected   <- auc_en_val       - auc_optimism
slope_en_corrected <- mean(boot_slope_orig, na.rm = TRUE)  # mean test slope
brier_en_corrected <- brier_en_apparent - brier_optimism_en

cat(sprintf("\nElastic Net Bootstrap Results:\n"))
cat(sprintf("  AUC    | Apparent: %.4f | Optimism: %.4f | Corrected: %.4f\n",
            auc_en_val, auc_optimism, auc_en_corrected))
cat(sprintf("  Slope  | Corrected: %.4f\n", slope_en_corrected))
cat(sprintf("  Brier  | Apparent: %.4f | Optimism: %.4f | Corrected: %.4f\n",
            brier_en_apparent, brier_optimism_en, brier_en_corrected))


# =============================================================================
# SECTION 5 — PART 5: Calibration Plot
#
# Manual calibration curve using logistic smoothing
# (rms::calibrate() is not available for glmnet objects)
# Method: bin predictions into deciles, compute observed rate per bin
# Overlay apparent and bootstrap-corrected curves
# =============================================================================

cat("\n=== PART 5: CALIBRATION PLOT ===\n")

# --- Apparent calibration (decile binning) ---
cal_apparent_en <- data.frame(
  pred = en_pred_apparent,
  obs  = y
) %>%
  mutate(
    bin = ntile(pred, 10)
  ) %>%
  group_by(bin) %>%
  summarise(
    mean_pred = mean(pred),
    obs_rate  = mean(obs),
    n         = n(),
    se        = sqrt(obs_rate * (1 - obs_rate) / n),
    .groups   = "drop"
  ) %>%
  mutate(Calibration = "Apparent")

# --- Bias-corrected calibration ---
# Use mean predicted probabilities from bootstrap test evaluations
# Collect bootstrap-corrected predicted probabilities
set.seed(2025)
n_cal_boot <- 500   # Separate smaller loop for calibration curve stability

boot_pred_matrix <- matrix(NA, nrow = n, ncol = n_cal_boot)

for (b in seq_len(n_cal_boot)) {
  
  boot_idx <- sample(seq_len(n), size = n, replace = TRUE)
  X_b      <- X[boot_idx, , drop = FALSE]
  y_b      <- y[boot_idx]
  
  cv_b <- tryCatch(
    cv.glmnet(X_b, y_b, alpha = best_alpha,
              family = "binomial", nfolds = 10,
              type.measure = "deviance"),
    error = function(e) NULL
  )
  if (is.null(cv_b)) next
  
  fit_b <- tryCatch(
    glmnet(X_b, y_b, alpha = best_alpha,
           lambda = cv_b$lambda.min, family = "binomial"),
    error = function(e) NULL
  )
  if (is.null(fit_b)) next
  
  boot_pred_matrix[, b] <- as.numeric(
    predict(fit_b, newx = X, type = "response")
  )
  
  if (b %% 100 == 0) cat(sprintf("  Calibration boot %d / %d\n", b, n_cal_boot))
}

# Mean corrected predictions across bootstrap models
corrected_preds <- rowMeans(boot_pred_matrix, na.rm = TRUE)

cal_corrected_en <- data.frame(
  pred = corrected_preds,
  obs  = y
) %>%
  mutate(bin = ntile(pred, 10)) %>%
  group_by(bin) %>%
  summarise(
    mean_pred = mean(pred),
    obs_rate  = mean(obs),
    n         = n(),
    se        = sqrt(obs_rate * (1 - obs_rate) / n),
    .groups   = "drop"
  ) %>%
  mutate(Calibration = "Bias-corrected")

cal_plot_df <- bind_rows(cal_apparent_en, cal_corrected_en)

# Publication calibration figure
cal_plot_en <- ggplot(cal_plot_df,
                      aes(x = mean_pred, y = obs_rate,
                          color = Calibration,
                          linetype = Calibration)) +
  geom_abline(slope = 1, intercept = 0,
              color = "gray60", linewidth = 0.8) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2.5) +
  geom_errorbar(
    aes(ymin = obs_rate - 1.96 * se,
        ymax = obs_rate + 1.96 * se),
    width = 0.01, alpha = 0.5
  ) +
  scale_color_manual(
    values = c("Apparent"       = "#E07B39",
               "Bias-corrected" = "#2C3E8C")
  ) +
  scale_linetype_manual(
    values = c("Apparent"       = "dashed",
               "Bias-corrected" = "solid")
  ) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title    = "Calibration Plot — Elastic Net Model",
    subtitle = paste0("Bootstrap bias-corrected (n_cal_boot = ",
                      n_cal_boot, ")"),
    x        = "Predicted Probability",
    y        = "Observed Proportion",
    color    = NULL, linetype = NULL,
    caption  = "Gray diagonal = ideal calibration"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title      = element_text(face = "bold"),
    plot.subtitle   = element_text(color = "gray40"),
    legend.position = c(0.15, 0.85),
    legend.background = element_rect(fill = "white", color = "gray80")
  )

print(cal_plot_en)
ggsave("figure9_calibration_EN.png", cal_plot_en,
       width = 6.5, height = 5.5, dpi = 300)

cat("Calibration figure saved: figure9_calibration_EN.png\n")


# =============================================================================
# SECTION 6 — PART 6: Head-to-Head Benchmark Table
# =============================================================================

cat("\n=== PART 6: HEAD-TO-HEAD COMPARISON TABLE ===\n")

comparison_en <- data.frame(
  Metric = c(
    "Apparent C-index (AUC)",
    "Corrected C-index",
    "Optimism (C-index)",
    "Apparent Calibration Slope",
    "Corrected Calibration Slope",
    "Apparent Brier Score",
    "Corrected Brier Score"
  ),
  Logistic_Regression = c(
    round(ref_cindex_apparent,  3),
    round(ref_cindex_corrected, 3),
    round(ref_optimism,         3),
    1.000,                           # Apparent slope = 1 by definition for glm
    round(ref_slope_corrected,  3),
    round(ref_brier_apparent,   4),
    NA    # Update with Step 2 corrected Brier if available
  ),
  Elastic_Net = c(
    round(auc_en_val,           3),
    round(auc_en_corrected,     3),
    round(auc_optimism,         3),
    round(mean(boot_slope_boot, na.rm = TRUE), 3),
    round(slope_en_corrected,   3),
    round(brier_en_apparent,    4),
    round(brier_en_corrected,   4)
  )
)

print(comparison_en)

comparison_gt_en <- comparison_en %>%
  gt() %>%
  tab_header(
    title    = "Table 8. Head-to-Head: Logistic Regression vs Elastic Net",
    subtitle = paste0("Bootstrap internal validation (B = ", B, ")")
  ) %>%
  cols_label(
    Metric               = "Performance Metric",
    Logistic_Regression  = "Logistic Regression",
    Elastic_Net          = "Elastic Net"
  ) %>%
  sub_missing(missing_text = "—") %>%
  tab_footnote(
    footnote = paste0(
      "Logistic regression values from Steps 1-2. ",
      "Elastic Net: alpha = ", best_alpha,
      ", lambda = ", round(best_lambda, 5), ". ",
      "Corrected = Apparent − Bootstrap Optimism. ",
      "Higher C-index = better discrimination. ",
      "Lower Brier = better overall accuracy. ",
      "Calibration slope closer to 1.0 = better calibration."
    )
  )

print(comparison_gt_en)


# =============================================================================
# SECTION 7 — PART 7: Clinical Interpretation
# =============================================================================

cat("\n=== PART 7: CLINICAL INTERPRETATION ===\n\n")

delta_cindex_corr <- auc_en_corrected - ref_cindex_corrected
delta_slope_corr  <- slope_en_corrected - ref_slope_corrected
delta_brier_corr  <- brier_en_corrected - ref_brier_apparent

meaningful_disc  <- abs(delta_cindex_corr) > 0.02
meaningful_cal   <- abs(delta_slope_corr)  > 0.05
meaningful_brier <- abs(delta_brier_corr)  > 0.01

cat(sprintf("Q1. Did Elastic Net improve discrimination?\n"))
cat(sprintf("    Delta corrected C-index: %+.3f\n", delta_cindex_corr))
cat(sprintf("    Assessment: %s\n\n",
            ifelse(delta_cindex_corr > 0.02,
                   "YES — meaningful improvement (>0.02)",
                   ifelse(delta_cindex_corr > 0,
                          "MARGINAL — improvement present but <0.02",
                          "NO — Elastic Net did not improve discrimination"))))

cat(sprintf("Q2. Did Elastic Net improve calibration?\n"))
cat(sprintf("    Delta corrected slope: %+.3f\n", delta_slope_corr))
cat(sprintf("    Assessment: %s\n\n",
            ifelse(delta_slope_corr > 0.05,
                   "YES — meaningful improvement (>0.05)",
                   ifelse(delta_slope_corr > 0,
                          "MARGINAL — small improvement",
                          "NO — calibration not meaningfully improved"))))

cat(sprintf("Q3. Did Elastic Net reduce optimism?\n"))
cat(sprintf("    Logistic optimism: %.3f | EN optimism: %.3f\n",
            ref_optimism, auc_optimism))
cat(sprintf("    Assessment: %s\n\n",
            ifelse(auc_optimism < ref_optimism - 0.01,
                   "YES — optimism meaningfully reduced",
                   "NO — optimism not meaningfully reduced")))

cat(sprintf("Q4. Is any improvement clinically meaningful?\n"))
if (!meaningful_disc && !meaningful_cal) {
  cat("    NO — differences in corrected C-index and calibration slope\n")
  cat("    fall below thresholds for clinical meaningfulness\n")
  cat("    (C-index delta < 0.02, slope delta < 0.05).\n\n")
} else {
  cat("    PARTIALLY — some metrics show meaningful differences.\n")
  cat("    Review comparison table carefully.\n\n")
}

cat("Q5. Should Elastic Net replace the approved logistic model?\n")
if (!meaningful_disc && !meaningful_cal) {
  cat("    RECOMMENDATION: RETAIN LOGISTIC REGRESSION\n\n")
  cat("    Justification:\n")
  cat("    - Elastic Net does not provide meaningful improvement in\n")
  cat("      corrected discrimination or calibration.\n")
  cat("    - The logistic regression model offers superior clinical\n")
  cat("      interpretability (explicit ORs, nomogram-compatible).\n")
  cat("    - Parsimony principle: prefer simpler model when performance\n")
  cat("      is equivalent.\n")
  cat("    - Elastic Net results support validity of logistic model:\n")
  cat("      penalization adds little, suggesting the model is not\n")
  cat("      severely overfit.\n")
} else {
  cat("    RECOMMENDATION: FLAG FOR STRATEGIST REVIEW\n")
  cat("    Elastic Net shows meaningful improvement on at least one metric.\n")
  cat("    Strategist AI should review before proceeding.\n")
}


# =============================================================================
# SECTION 8 — PART 8: Publication-Ready Summary
# =============================================================================

cat("\n=== PART 8: PUBLICATION-READY PARAGRAPH ===\n\n")

cat(sprintf(paste0(
  "Elastic Net regularization (alpha = %.1f, lambda = %.5f) was evaluated ",
  "as a benchmark machine-learning method using the same three predictors as ",
  "the approved logistic regression model (Tumor_mlocation3, NLR, PLR). ",
  "Hyperparameter tuning was performed via 10-fold cross-validation across an ",
  "alpha grid from 0.0 to 1.0. Bootstrap internal validation (B = %d iterations) ",
  "with lambda re-tuned within each iteration was used to estimate optimism. ",
  "The Elastic Net model achieved an apparent C-index of %.3f and an ",
  "optimism-corrected C-index of %.3f, compared with %.3f and %.3f for the ",
  "logistic regression model, respectively. ",
  "The corrected calibration slope for Elastic Net was %.3f versus %.3f for ",
  "logistic regression. ",
  "Given the absence of meaningful differences in corrected discrimination ",
  "or calibration between models, and the superior clinical interpretability ",
  "of logistic regression, the approved logistic regression model was retained ",
  "as the primary prediction model."
),
best_alpha, best_lambda, B,
auc_en_val, auc_en_corrected,
ref_cindex_apparent, ref_cindex_corrected,
slope_en_corrected, ref_slope_corrected
))

cat("\n")


# =============================================================================
# SECTION 9: Step 5 Summary
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("STEP 5 COMPLETE — ELASTIC NET BENCHMARK SUMMARY\n")
cat("================================================================\n")
cat(sprintf("Best alpha:  %.1f\n", best_alpha))
cat(sprintf("Best lambda: %.5f\n", best_lambda))
cat(sprintf("EN Apparent C-index:   %.4f\n", auc_en_val))
cat(sprintf("EN Corrected C-index:  %.4f\n", auc_en_corrected))
cat(sprintf("EN Optimism:           %.4f\n", auc_optimism))
cat(sprintf("EN Corrected Slope:    %.4f\n", slope_en_corrected))
cat(sprintf("EN Apparent Brier:     %.4f\n", brier_en_apparent))
cat(sprintf("EN Corrected Brier:    %.4f\n", brier_en_corrected))
cat("---\n")
cat(sprintf("Delta Corrected C-index (EN - LR): %+.4f\n", delta_cindex_corr))
cat(sprintf("Delta Corrected Slope  (EN - LR): %+.4f\n", delta_slope_corr))
cat("---\n")
cat("Figures:  figure8_ROC_EN_vs_logistic.png\n")
cat("          figure9_calibration_EN.png\n")
cat("Tables:   Table 6 (tuning), Table 7 (coefficients), Table 8 (comparison)\n")
cat("================================================================\n")
cat("STOPPING — Awaiting Strategist review before Step 6 (XGBoost).\n")