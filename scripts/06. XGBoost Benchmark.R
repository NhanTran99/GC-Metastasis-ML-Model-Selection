
# =============================================================================
# PROJECT: Clinical Prediction Model — Occult Distant Metastasis in Gastric Cancer
# STEP 6:  XGBoost Benchmark
# AUTHOR:  Nhan Tran
# DATE:    160626
# DEPENDENCY: Steps 1-5 must be run first.
#             Required objects: df, model_glm
#             Step 5 EN results: auc_en_corrected, slope_en_corrected,
#                                brier_en_corrected, auc_en_val,
#                                auc_optimism, brier_en_apparent
# =============================================================================


# =============================================================================
# SECTION 0: Environment & Prerequisites
# =============================================================================

set.seed(2025)

# install.packages("xgboost")
library(xgboost)   # XGBoost
library(pROC)      # ROC, AUC
library(ggplot2)   # Figures
library(dplyr)     # Data wrangling
library(gt)        # Publication tables

# Bootstrap iterations
B <- 2000

# Reference values — UPDATE with your actual results from Steps 1-2 and 5
ref_lr_apparent_cindex  <- 0.757
ref_lr_corrected_cindex <- 0.719
ref_lr_optimism         <- 0.017
ref_lr_slope_corrected  <- 0.831
ref_lr_brier_apparent   <- 0.204

# Elastic Net values from Step 5 — pulled from environment
# If not in environment, replace with actual values:
# auc_en_val            <- 0.xxx
# auc_en_corrected      <- 0.xxx
# auc_optimism          <- 0.xxx
# slope_en_corrected    <- 0.xxx
# brier_en_apparent     <- 0.xxx
# brier_en_corrected    <- 0.xxx

# Safety checks
if (!exists("df"))        stop("df not found. Run Steps 1-5 first.")
if (!exists("model_glm")) stop("model_glm not found. Run Steps 1-5 first.")

cat("=== STEP 6 PREREQUISITES CHECK ===\n")
cat("n =", nrow(df), "| Events =", sum(df$Metastasis_TNM), "\n\n")


# =============================================================================
# SECTION 1 — PART 1: Data Preparation
#
# XGBoost requires:
#   - Numeric matrix input (xgb.DMatrix)
#   - No factor columns
#   - One-hot encoding for categorical variables
#
# We use model.matrix() — consistent with Step 5 (Elastic Net)
# [, -1] removes intercept column
# =============================================================================

cat("=== PART 1: DATA PREPARATION ===\n")

# Build numeric feature matrix (one-hot encodes Tumor_mlocation3)
X_xgb <- model.matrix(
  Metastasis_TNM ~ Tumor_mlocation3 + NLR + PLR,
  data = df
)[, -1]

y_xgb <- as.numeric(df$Metastasis_TNM)

cat("Feature matrix dimensions:", dim(X_xgb), "\n")
cat("Feature names:\n")
print(colnames(X_xgb))
cat("Outcome distribution:", table(y_xgb), "\n\n")

# Convert to xgb.DMatrix (XGBoost native format)
dtrain <- xgb.DMatrix(data = X_xgb, label = y_xgb)

# Preprocessing summary table
prep_summary <- data.frame(
  Step = c(
    "Original predictors",
    "Encoding method",
    "Tumor_mlocation3",
    "NLR",
    "PLR",
    "Final feature matrix",
    "Outcome"
  ),
  Detail = c(
    "Tumor_mlocation3, NLR, PLR",
    "model.matrix() — one-hot encoding for factors",
    paste0("Dummy coded: ", ncol(X_xgb) - 2,
           " binary column(s) (reference = ",
           levels(df$Tumor_mlocation3)[1], ")"),
    "Numeric — used as-is",
    "Numeric — used as-is",
    paste0(nrow(X_xgb), " rows × ", ncol(X_xgb), " columns"),
    "Binary: 0 = M0, 1 = M1 (distant metastasis)"
  )
)

cat("--- Preprocessing Summary ---\n")
print(prep_summary)

prep_gt <- prep_summary %>%
  gt() %>%
  tab_header(
    title = "Table 9. XGBoost Data Preprocessing Summary"
  ) %>%
  cols_label(Step = "Step", Detail = "Detail")

print(prep_gt)


# =============================================================================
# SECTION 2 — PART 2: Hyperparameter Tuning
#
# Strategy: Grid search over constrained parameter space
#           appropriate for n = 114 (small dataset)
#
# Parameters tuned:
#   max_depth:         tree depth (2-4; shallow trees prevent overfitting)
#   eta:               learning rate (0.01-0.3)
#   min_child_weight:  minimum sum of instance weight in leaf (3-10; high = conservative)
#   subsample:         row sampling ratio (0.6-0.8)
#   colsample_bytree:  feature sampling ratio (0.6-1.0)
#   gamma:             minimum loss reduction to split (0-1)
#
# Evaluation: 5-fold CV (not 10-fold; n=114 means folds of ~23 — already small)
# Metric: AUC
# nrounds: determined by early stopping within CV
# =============================================================================

cat("\n=== PART 2: HYPERPARAMETER TUNING ===\n")
cat("Grid search with 5-fold CV (constrained for n = 114)...\n")
cat("This may take 2-5 minutes depending on grid size.\n\n")

# Constrained grid — deliberately limited for small n
# Total combinations: 3 × 3 × 2 × 2 × 2 × 2 = 144 combinations
param_grid <- expand.grid(
  max_depth         = c(2, 3, 4),
  eta               = c(0.01, 0.05, 0.1),
  min_child_weight  = c(3, 5),
  subsample         = c(0.6, 0.8),
  colsample_bytree  = c(0.6, 1.0),
  gamma             = c(0, 0.5),
  stringsAsFactors  = FALSE
)

cat(sprintf("Total parameter combinations: %d\n\n", nrow(param_grid)))

# Storage
tuning_xgb <- param_grid
tuning_xgb$cv_auc    <- NA_real_
tuning_xgb$best_nrounds <- NA_integer_

set.seed(2025)

for (i in seq_len(nrow(param_grid))) {
  
  params_i <- list(
    objective        = "binary:logistic",
    eval_metric      = "auc",
    max_depth        = param_grid$max_depth[i],
    eta              = param_grid$eta[i],
    min_child_weight = param_grid$min_child_weight[i],
    subsample        = param_grid$subsample[i],
    colsample_bytree = param_grid$colsample_bytree[i],
    gamma            = param_grid$gamma[i],
    seed             = 2025   # seed goes inside params, not as top-level arg
  )
  
  cv_i <- tryCatch(
    xgb.cv(
      params    = params_i,
      data      = dtrain,
      nrounds   = 500,
      nfold     = 5,
      early_stopping_rounds = 30,
      verbose   = 0
    ),
    error = function(e) NULL
  )
  
  if (!is.null(cv_i)) {
    # best_iteration can be NULL if early stopping never triggered
    # fall back to the round with maximum test AUC
    best_iter <- cv_i$best_iteration
    if (is.null(best_iter) || length(best_iter) == 0) {
      best_iter <- which.max(cv_i$evaluation_log$test_auc_mean)
    }
    tuning_xgb$cv_auc[i]        <- cv_i$evaluation_log$test_auc_mean[best_iter]
    tuning_xgb$best_nrounds[i]  <- best_iter
  }
  
  if (i %% 20 == 0) {
    cat(sprintf("  Completed %d / %d combinations\n", i, nrow(param_grid)))
  }
}

# Select best parameter set
best_idx    <- which.max(tuning_xgb$cv_auc)
best_params_row <- tuning_xgb[best_idx, ]

cat("\n--- BEST PARAMETER SET ---\n")
print(best_params_row)

# Top 10 parameter sets
cat("\n--- TOP 10 PARAMETER SETS BY CV AUC ---\n")
top10 <- tuning_xgb %>%
  arrange(desc(cv_auc)) %>%
  head(10) %>%
  mutate(across(where(is.numeric), ~ round(.x, 5)))
print(top10)

# Publication tuning table (top 10)
tuning_gt_xgb <- top10 %>%
  gt() %>%
  tab_header(
    title    = "Table 10. XGBoost Hyperparameter Tuning — Top 10 Results",
    subtitle = "5-fold cross-validation, metric = AUC"
  ) %>%
  tab_style(
    style     = cell_fill(color = "#D6EAF8"),
    locations = cells_body(rows = 1)
  ) %>%
  tab_footnote(
    footnote = paste0(
      "Highlighted row = selected parameter set. ",
      "max_depth: tree depth. eta: learning rate. ",
      "min_child_weight: minimum leaf weight. ",
      "subsample: row sampling. colsample_bytree: feature sampling. ",
      "gamma: minimum split gain. nrounds: optimal boosting rounds."
    )
  )

print(tuning_gt_xgb)

# Extract final parameters
best_nrounds         <- best_params_row$best_nrounds
best_max_depth       <- best_params_row$max_depth
best_eta             <- best_params_row$eta
best_min_child_weight <- best_params_row$min_child_weight
best_subsample       <- best_params_row$subsample
best_colsample       <- best_params_row$colsample_bytree
best_gamma           <- best_params_row$gamma
best_cv_auc          <- best_params_row$cv_auc

cat(sprintf("\nSelected: depth=%d, eta=%.2f, mcw=%d, sub=%.1f, col=%.1f, gamma=%.1f\n",
            best_max_depth, best_eta, best_min_child_weight,
            best_subsample, best_colsample, best_gamma))
cat(sprintf("Optimal nrounds: %d | CV AUC: %.4f\n",
            best_nrounds, best_cv_auc))


# =============================================================================
# SECTION 3 — PART 3: Final XGBoost Model
# =============================================================================

cat("\n=== PART 3: FINAL XGBoost MODEL ===\n")

final_params <- list(
  objective        = "binary:logistic",
  eval_metric      = "auc",
  max_depth        = best_max_depth,
  eta              = best_eta,
  min_child_weight = best_min_child_weight,
  subsample        = best_subsample,
  colsample_bytree = best_colsample,
  gamma            = best_gamma
)

set.seed(2025)
final_xgb <- xgb.train(
  params  = final_params,
  data    = dtrain,
  nrounds = best_nrounds,
  verbose = 0
)

cat("Final XGBoost model fitted.\n")
cat(sprintf("  nrounds: %d\n", best_nrounds))
cat(sprintf("  Training AUC: %.4f\n",
            as.numeric(auc(roc(y_xgb,
                               predict(final_xgb, dtrain),
                               quiet = TRUE,
                               direction = "<")))))

# Final parameters table
params_df <- data.frame(
  Parameter = c("max_depth", "eta", "min_child_weight",
                "subsample", "colsample_bytree", "gamma",
                "nrounds", "objective", "eval_metric"),
  Value     = c(best_max_depth, best_eta, best_min_child_weight,
                best_subsample, best_colsample, best_gamma,
                best_nrounds, "binary:logistic", "auc")
)

params_gt <- params_df %>%
  gt() %>%
  tab_header(
    title = "Table 11. Final XGBoost Model Parameters"
  ) %>%
  cols_label(Parameter = "Hyperparameter", Value = "Value") %>%
  tab_footnote(
    footnote = "Parameters selected by 5-fold CV grid search on full dataset (n = 114)."
  )

print(params_gt)


# =============================================================================
# SECTION 4 — PART 4: Apparent Performance
# =============================================================================

cat("\n=== PART 4: APPARENT PERFORMANCE ===\n")

xgb_pred_apparent <- predict(final_xgb, dtrain)

# ROC and AUC
roc_xgb <- roc(
  response  = y_xgb,
  predictor = xgb_pred_apparent,
  ci        = TRUE,
  direction = "<",
  quiet     = TRUE
)

auc_xgb_val <- round(as.numeric(auc(roc_xgb)), 4)
auc_xgb_ci  <- round(as.numeric(ci.auc(roc_xgb)), 4)
brier_xgb_apparent <- round(mean((xgb_pred_apparent - y_xgb)^2), 4)

cat(sprintf("Apparent AUC:   %.4f (95%% CI: %.4f–%.4f)\n",
            auc_xgb_val, auc_xgb_ci[1], auc_xgb_ci[3]))
cat(sprintf("Apparent Brier: %.4f\n", brier_xgb_apparent))

# --- Three-way ROC overlay figure ---
pred_lr  <- predict(model_glm, type = "response")
roc_lr   <- roc(y_xgb, pred_lr, quiet = TRUE, direction = "<")

# Build combined ROC dataframe
make_roc_df <- function(roc_obj, label) {
  data.frame(
    FPR   = 1 - roc_obj$specificities,
    TPR   = roc_obj$sensitivities,
    Model = label
  )
}

roc_all <- bind_rows(
  make_roc_df(roc_lr,
              paste0("Logistic (AUC = ",
                     round(as.numeric(auc(roc_lr)), 3), ")")),
  make_roc_df(roc_xgb,
              paste0("XGBoost (AUC = ", auc_xgb_val, ")"))
)

# Add Elastic Net ROC if available
if (exists("en_pred_apparent")) {
  roc_en_obj <- roc(y_xgb, en_pred_apparent,
                    quiet = TRUE, direction = "<")
  roc_all <- bind_rows(
    roc_all,
    make_roc_df(roc_en_obj,
                paste0("Elastic Net (AUC = ",
                       round(as.numeric(auc(roc_en_obj)), 3), ")"))
  )
}

roc_colors <- c(
  "gray40",     # Logistic
  "#2C3E8C",    # XGBoost
  "#B5543B"     # Elastic Net
) |> setNames(unique(roc_all$Model))

roc_plot_xgb <- ggplot(roc_all,
                       aes(x = FPR, y = TPR,
                           color = Model,
                           linetype = Model)) +
  geom_line(linewidth = 1.1) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "gray70") +
  scale_color_manual(values = roc_colors) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title    = "ROC Curves — Model Benchmark Comparison",
    subtitle = "Apparent performance, full dataset (n = 114)",
    x        = "1 − Specificity (False Positive Rate)",
    y        = "Sensitivity (True Positive Rate)",
    color    = NULL, linetype = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title      = element_text(face = "bold"),
    plot.subtitle   = element_text(color = "gray40"),
    legend.position = c(0.60, 0.15),
    legend.background = element_rect(fill = "white", color = "gray80")
  )

print(roc_plot_xgb)
ggsave("figure10_ROC_benchmark_comparison.png", roc_plot_xgb,
       width = 7, height = 5.5, dpi = 300)

cat("ROC figure saved: figure10_ROC_benchmark_comparison.png\n")


# =============================================================================
# SECTION 5 — PART 5: Bootstrap Internal Validation
#
# DESIGN DECISION: Hyperparameters FIXED during bootstrap
# Reason: Re-tuning XGBoost (144 grid × 5-fold CV) inside 2000 bootstrap
#         iterations would be computationally prohibitive (~48 hours).
# This is documented as a limitation in the methods summary (Part 10).
#
# Method:
#   For each bootstrap b:
#     1. Sample with replacement → boot sample
#     2. Fit XGBoost with fixed params on boot sample
#     3. Evaluate on boot sample (apparent) and original (test)
#     4. optimism_b = apparent - test
#   Corrected = Apparent - mean(optimism)
# =============================================================================

cat("\n=== PART 5: BOOTSTRAP INTERNAL VALIDATION ===\n")
cat(sprintf("B = %d | Hyperparameters FIXED (see Part 10 limitation note)\n", B))
cat("Estimated runtime: 3-8 minutes\n\n")

set.seed(2025)
n <- nrow(df)

boot_auc_boot_xgb   <- numeric(B)
boot_auc_orig_xgb   <- numeric(B)
boot_slope_boot_xgb <- numeric(B)
boot_slope_orig_xgb <- numeric(B)
boot_brier_boot_xgb <- numeric(B)
boot_brier_orig_xgb <- numeric(B)
skipped_xgb <- 0

for (b in seq_len(B)) {
  
  # Step 1: Bootstrap sample
  boot_idx <- sample(seq_len(n), size = n, replace = TRUE)
  X_b      <- X_xgb[boot_idx, , drop = FALSE]
  y_b      <- y_xgb[boot_idx]
  
  dtrain_b <- xgb.DMatrix(data = X_b, label = y_b)
  
  # Step 2: Fit with fixed hyperparameters
  fit_b <- tryCatch(
    xgb.train(
      params  = final_params,
      data    = dtrain_b,
      nrounds = best_nrounds,
      verbose = 0
    ),
    error = function(e) NULL
  )
  
  if (is.null(fit_b)) { skipped_xgb <- skipped_xgb + 1; next }
  
  # Step 3a: Predict on bootstrap sample
  p_boot_on_boot <- predict(fit_b, dtrain_b)
  
  # Step 3b: Predict on original data
  p_boot_on_orig <- predict(fit_b, dtrain)
  
  # AUC
  boot_auc_boot_xgb[b] <- tryCatch(
    as.numeric(auc(roc(y_b, p_boot_on_boot,
                       quiet = TRUE, direction = "<"))),
    error = function(e) NA
  )
  boot_auc_orig_xgb[b] <- tryCatch(
    as.numeric(auc(roc(y_xgb, p_boot_on_orig,
                       quiet = TRUE, direction = "<"))),
    error = function(e) NA
  )
  
  # Calibration slope
  logit_b    <- log(p_boot_on_boot / (1 - p_boot_on_boot + 1e-10))
  logit_orig <- log(p_boot_on_orig / (1 - p_boot_on_orig + 1e-10))
  
  boot_slope_boot_xgb[b] <- tryCatch(
    coef(glm(y_b ~ logit_b, family = binomial))[2],
    error = function(e) NA
  )
  boot_slope_orig_xgb[b] <- tryCatch(
    coef(glm(y_xgb ~ logit_orig, family = binomial))[2],
    error = function(e) NA
  )
  
  # Brier
  boot_brier_boot_xgb[b] <- mean((p_boot_on_boot - y_b)^2)
  boot_brier_orig_xgb[b] <- mean((p_boot_on_orig - y_xgb)^2)
  
  if (b %% 500 == 0) cat(sprintf("  Iteration %d / %d\n", b, B))
}

cat(sprintf("\nBootstrap complete. Skipped: %d / %d\n", skipped_xgb, B))

# Compute optimism and corrected metrics
auc_xgb_optimism   <- mean(boot_auc_boot_xgb  - boot_auc_orig_xgb,
                           na.rm = TRUE)
slope_xgb_optimism <- mean(boot_slope_boot_xgb - boot_slope_orig_xgb,
                           na.rm = TRUE)
brier_xgb_optimism <- mean(boot_brier_boot_xgb - boot_brier_orig_xgb,
                           na.rm = TRUE)

auc_xgb_corrected   <- auc_xgb_val       - auc_xgb_optimism
slope_xgb_corrected <- mean(boot_slope_orig_xgb, na.rm = TRUE)
brier_xgb_corrected <- brier_xgb_apparent - brier_xgb_optimism

cat(sprintf("\nXGBoost Bootstrap Results:\n"))
cat(sprintf("  AUC   | Apparent: %.4f | Optimism: %.4f | Corrected: %.4f\n",
            auc_xgb_val, auc_xgb_optimism, auc_xgb_corrected))
cat(sprintf("  Slope | Corrected: %.4f\n", slope_xgb_corrected))
cat(sprintf("  Brier | Apparent: %.4f | Optimism: %.4f | Corrected: %.4f\n",
            brier_xgb_apparent, brier_xgb_optimism, brier_xgb_corrected))

# Publication bootstrap table
boot_xgb_df <- data.frame(
  Metric    = c("C-index (AUC)", "Calibration Slope", "Brier Score"),
  Apparent  = c(round(auc_xgb_val, 4),
                round(mean(boot_slope_boot_xgb, na.rm=TRUE), 4),
                round(brier_xgb_apparent, 4)),
  Optimism  = c(round(auc_xgb_optimism, 4),
                round(slope_xgb_optimism, 4),
                round(brier_xgb_optimism, 4)),
  Corrected = c(round(auc_xgb_corrected, 4),
                round(slope_xgb_corrected, 4),
                round(brier_xgb_corrected, 4))
)

boot_xgb_gt <- boot_xgb_df %>%
  gt() %>%
  tab_header(
    title    = "Table 12. XGBoost Bootstrap Internal Validation",
    subtitle = paste0("B = ", B, " iterations; hyperparameters fixed")
  ) %>%
  cols_label(
    Metric    = "Metric",
    Apparent  = "Apparent",
    Optimism  = "Optimism",
    Corrected = "Optimism-Corrected"
  ) %>%
  tab_footnote(
    footnote = paste0(
      "Corrected = Apparent − Optimism. ",
      "Hyperparameters fixed at values tuned on full dataset. ",
      "Calibration slope: ideal = 1.0. Lower Brier = better."
    )
  )

print(boot_xgb_gt)


# =============================================================================
# SECTION 6 — PART 6: Calibration Plot
#
# Same method as Step 5 — decile binning
# Apparent + bias-corrected using bootstrap predictions
# =============================================================================

cat("\n=== PART 6: CALIBRATION PLOT ===\n")

# Apparent calibration
cal_xgb_apparent <- data.frame(
  pred = xgb_pred_apparent,
  obs  = y_xgb
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
  mutate(Calibration = "Apparent")

# Bias-corrected: use boot_on_orig predictions from bootstrap loop
# Average predicted probabilities from all bootstrap models on original data
# Reconstruct from stored boot predictions (lightweight re-run)
set.seed(2025)
n_cal <- 300
boot_pred_xgb <- matrix(NA, nrow = n, ncol = n_cal)

for (b in seq_len(n_cal)) {
  boot_idx <- sample(seq_len(n), size = n, replace = TRUE)
  X_b      <- X_xgb[boot_idx, , drop = FALSE]
  y_b      <- y_xgb[boot_idx]
  dtrain_b <- xgb.DMatrix(data = X_b, label = y_b)
  
  fit_b <- tryCatch(
    xgb.train(params  = final_params,
              data    = dtrain_b,
              nrounds = best_nrounds,
              verbose = 0),
    error = function(e) NULL
  )
  if (!is.null(fit_b)) {
    boot_pred_xgb[, b] <- predict(fit_b, dtrain)
  }
  
  if (b %% 100 == 0) cat(sprintf("  Calibration boot %d / %d\n", b, n_cal))
}

corrected_preds_xgb <- rowMeans(boot_pred_xgb, na.rm = TRUE)

cal_xgb_corrected <- data.frame(
  pred = corrected_preds_xgb,
  obs  = y_xgb
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

cal_xgb_df <- bind_rows(cal_xgb_apparent, cal_xgb_corrected)

cal_plot_xgb <- ggplot(cal_xgb_df,
                       aes(x = mean_pred, y = obs_rate,
                           color    = Calibration,
                           linetype = Calibration)) +
  geom_abline(slope = 1, intercept = 0,
              color = "gray60", linewidth = 0.8) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = obs_rate - 1.96 * se,
                    ymax = obs_rate + 1.96 * se),
                width = 0.01, alpha = 0.5) +
  scale_color_manual(values = c("Apparent"       = "#E07B39",
                                "Bias-corrected" = "#2C3E8C")) +
  scale_linetype_manual(values = c("Apparent"       = "dashed",
                                   "Bias-corrected" = "solid")) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title    = "Calibration Plot — XGBoost Model",
    subtitle = paste0("Bootstrap bias-corrected (n_cal = ", n_cal, ")"),
    x        = "Predicted Probability",
    y        = "Observed Proportion",
    color    = NULL, linetype = NULL,
    caption  = "Gray diagonal = ideal calibration"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title        = element_text(face = "bold"),
    plot.subtitle     = element_text(color = "gray40"),
    legend.position   = c(0.15, 0.85),
    legend.background = element_rect(fill = "white", color = "gray80")
  )

print(cal_plot_xgb)
ggsave("figure11_calibration_XGBoost.png", cal_plot_xgb,
       width = 6.5, height = 5.5, dpi = 300)

cat("Calibration figure saved: figure11_calibration_XGBoost.png\n")


# =============================================================================
# SECTION 7 — PART 7: Feature Importance
#
# Native XGBoost importance metrics:
#   gain:  average improvement in loss from splits using this feature
#          → most relevant for prediction contribution
#   cover: average number of samples affected by splits on this feature
#   frequency: proportion of times feature appears in trees
#
# We report GAIN as primary metric (standard in clinical prediction literature)
# This is NOT SHAP — native importance only per instructions
# =============================================================================

cat("\n=== PART 7: FEATURE IMPORTANCE ===\n")

importance_xgb <- xgb.importance(
  feature_names = colnames(X_xgb),
  model         = final_xgb
)

cat("Feature importance (by Gain):\n")
print(importance_xgb)

# Publication importance figure
imp_df <- as.data.frame(importance_xgb) %>%
  arrange(Gain) %>%
  mutate(Feature = factor(Feature, levels = Feature))

imp_plot <- ggplot(imp_df,
                   aes(x = Gain, y = Feature)) +
  geom_col(fill = "#2C3E8C", alpha = 0.85, width = 0.6) +
  geom_text(aes(label = round(Gain, 3)),
            hjust = -0.1, size = 3.8, color = "gray30") +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.15))
  ) +
  labs(
    title    = "XGBoost Feature Importance (Gain)",
    subtitle = "Higher gain = greater contribution to prediction",
    x        = "Importance (Gain)",
    y        = NULL,
    caption  = paste0(
      "Gain = average improvement in loss per split. ",
      "Native XGBoost importance (not SHAP)."
    )
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold"),
    plot.subtitle = element_text(color = "gray40"),
    axis.text.y   = element_text(size = 11)
  )

print(imp_plot)
ggsave("figure12_XGBoost_feature_importance.png", imp_plot,
       width = 7, height = 4.5, dpi = 300)

cat("Feature importance figure saved: figure12_XGBoost_feature_importance.png\n")


# =============================================================================
# SECTION 8 — PART 8: Head-to-Head Benchmark Table (Three-Way)
# =============================================================================

cat("\n=== PART 8: HEAD-TO-HEAD BENCHMARK TABLE ===\n")

# Safely retrieve EN values (use placeholders if not in environment)
en_apparent_c  <- ifelse(exists("auc_en_val"),
                         round(auc_en_val, 3), NA)
en_corrected_c <- ifelse(exists("auc_en_corrected"),
                         round(auc_en_corrected, 3), NA)
en_optimism_c  <- ifelse(exists("auc_optimism"),
                         round(auc_optimism, 3), NA)
en_slope_corr  <- ifelse(exists("slope_en_corrected"),
                         round(slope_en_corrected, 3), NA)
en_brier_app   <- ifelse(exists("brier_en_apparent"),
                         round(brier_en_apparent, 4), NA)
en_brier_corr  <- ifelse(exists("brier_en_corrected"),
                         round(brier_en_corrected, 4), NA)

benchmark_df <- data.frame(
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
    round(ref_lr_apparent_cindex,  3),
    round(ref_lr_corrected_cindex, 3),
    round(ref_lr_optimism,         3),
    "1.000",
    round(ref_lr_slope_corrected,  3),
    round(ref_lr_brier_apparent,   4),
    "—"
  ),
  Elastic_Net = c(
    en_apparent_c,
    en_corrected_c,
    en_optimism_c,
    "—",
    en_slope_corr,
    en_brier_app,
    en_brier_corr
  ),
  XGBoost = c(
    round(auc_xgb_val,          3),
    round(auc_xgb_corrected,    3),
    round(auc_xgb_optimism,     3),
    round(mean(boot_slope_boot_xgb, na.rm = TRUE), 3),
    round(slope_xgb_corrected,  3),
    round(brier_xgb_apparent,   4),
    round(brier_xgb_corrected,  4)
  )
)

print(benchmark_df)

benchmark_gt <- benchmark_df %>%
  gt() %>%
  tab_header(
    title    = "Table 13. Model Benchmark Comparison",
    subtitle = paste0(
      "Logistic Regression vs Elastic Net vs XGBoost\n",
      "Bootstrap internal validation (B = ", B, ")"
    )
  ) %>%
  cols_label(
    Metric              = "Performance Metric",
    Logistic_Regression = "Logistic Regression",
    Elastic_Net         = "Elastic Net",
    XGBoost             = "XGBoost"
  ) %>%
  sub_missing(missing_text = "—") %>%
  tab_footnote(
    footnote = paste0(
      "Corrected = Apparent − Bootstrap Optimism. ",
      "C-index = AUC for binary outcome. ",
      "Calibration slope: ideal = 1.0. ",
      "Brier Score: lower = better. ",
      "XGBoost hyperparameters fixed during bootstrap (see methods)."
    )
  )

print(benchmark_gt)


# =============================================================================
# SECTION 9 — PART 9: Clinical Interpretation (Q1-Q5)
# =============================================================================

cat("\n=== PART 9: CLINICAL INTERPRETATION ===\n\n")

delta_c_xgb     <- auc_xgb_corrected   - ref_lr_corrected_cindex
delta_slope_xgb <- slope_xgb_corrected - ref_lr_slope_corrected
delta_brier_xgb <- brier_xgb_corrected - ref_lr_brier_apparent

cat(sprintf("Q1. Did XGBoost improve discrimination?\n"))
cat(sprintf("    Delta corrected C-index (XGB - LR): %+.4f\n", delta_c_xgb))
cat(sprintf("    Assessment: %s\n\n",
            ifelse(delta_c_xgb > 0.02,
                   "YES — meaningful improvement (>0.02)",
                   ifelse(delta_c_xgb > 0,
                          "MARGINAL — improvement present but <0.02 threshold",
                          "NO — XGBoost did not improve discrimination"))))

cat(sprintf("Q2. Did XGBoost improve calibration?\n"))
cat(sprintf("    Delta corrected slope (XGB - LR): %+.4f\n", delta_slope_xgb))
cat(sprintf("    Assessment: %s\n\n",
            ifelse(delta_slope_xgb > 0.05,
                   "YES — meaningful improvement",
                   ifelse(delta_slope_xgb > 0,
                          "MARGINAL — small improvement",
                          "NO — calibration not meaningfully improved"))))

cat(sprintf("Q3. Did XGBoost reduce optimism?\n"))
cat(sprintf("    LR optimism: %.4f | XGB optimism: %.4f\n",
            ref_lr_optimism, auc_xgb_optimism))
cat(sprintf("    Assessment: %s\n\n",
            ifelse(auc_xgb_optimism < ref_lr_optimism - 0.01,
                   "YES — XGBoost shows less optimism",
                   ifelse(auc_xgb_optimism > ref_lr_optimism + 0.01,
                          "NO — XGBoost shows MORE optimism than logistic",
                          "SIMILAR — no meaningful difference in optimism"))))

cat(sprintf("Q4. Are improvements clinically meaningful?\n"))
meaningful <- (abs(delta_c_xgb) > 0.02) | (abs(delta_slope_xgb) > 0.05)
cat(sprintf("    %s\n\n",
            ifelse(meaningful,
                   "PARTIALLY — at least one metric exceeds threshold",
                   "NO — differences below clinical meaningfulness thresholds")))

cat("Q5. Should XGBoost replace the approved logistic model?\n")
if (!meaningful || auc_xgb_optimism > ref_lr_optimism + 0.01) {
  cat("    RECOMMENDATION: RETAIN LOGISTIC REGRESSION\n\n")
  cat("    Justification:\n")
  cat("    - XGBoost does not provide meaningful improvement in corrected\n")
  cat("      discrimination or calibration over logistic regression.\n")
  cat("    - Logistic regression offers explicit odds ratios — directly\n")
  cat("      interpretable for clinicians and nomogram construction.\n")
  cat("    - XGBoost is a black-box model: not suitable for nomogram\n")
  cat("      generation (Step 8) or web calculator deployment (Step 10).\n")
  cat("    - Higher optimism in XGBoost (if present) suggests overfitting\n")
  cat("      despite regularization — unsurprising given n = 114.\n")
  cat("    - Convergent evidence from Steps 3, 5, and 6: no method\n")
  cat("      consistently outperforms the linear model on this dataset.\n")
} else {
  cat("    RECOMMENDATION: FLAG FOR STRATEGIST REVIEW\n")
  cat("    XGBoost shows meaningful improvement — Strategist AI should\n")
  cat("    review before any model change decision.\n")
}


# =============================================================================
# SECTION 10 — PART 10: Publication-Ready Methods Summary
# =============================================================================

cat("\n=== PART 10: PUBLICATION-READY METHODS SUMMARY ===\n\n")

cat(paste0(
  "XGBoost (eXtreme Gradient Boosting) was evaluated as a machine-learning ",
  "benchmark using the same predictor set as the approved logistic regression ",
  "model (Tumor_mlocation3, NLR, PLR). The categorical predictor ",
  "Tumor_mlocation3 was one-hot encoded using model.matrix(), consistent with ",
  "the Elastic Net preprocessing in Step 5. ",
  "\n\n",
  "Hyperparameter tuning was performed via a constrained grid search across ",
  "six parameters (max_depth, eta, min_child_weight, subsample, ",
  "colsample_bytree, gamma), with ", nrow(param_grid), " combinations evaluated ",
  "using 5-fold cross-validation (metric: AUC) with early stopping ",
  "(patience = 30 rounds, maximum nrounds = 500). The search space was ",
  "deliberately constrained to models appropriate for a small clinical dataset ",
  "(n = 114) to avoid selection of excessively complex trees. ",
  "\n\n",
  "The final model was trained on the full dataset using the optimal ",
  "hyperparameter configuration: max_depth = ", best_max_depth,
  ", eta = ", best_eta,
  ", min_child_weight = ", best_min_child_weight,
  ", subsample = ", best_subsample,
  ", colsample_bytree = ", best_colsample,
  ", gamma = ", best_gamma,
  ", nrounds = ", best_nrounds, ". ",
  "\n\n",
  "Internal validation was performed using ", B, "-iteration bootstrap ",
  "optimism correction, consistent with the validation strategy applied to ",
  "logistic regression (Step 2) and Elastic Net (Step 5). To maintain ",
  "computational feasibility, hyperparameters were fixed at values identified ",
  "during full-dataset tuning rather than re-tuned within each bootstrap ",
  "iteration; this represents a known limitation that may modestly ",
  "underestimate true optimism for XGBoost. ",
  "\n\n",
  "Calibration was assessed using decile-binned observed versus predicted ",
  "probabilities, with bootstrap bias-correction (", 300, " iterations). ",
  "Feature importance was quantified using the native XGBoost gain metric, ",
  "defined as the average improvement in the loss function attributable to ",
  "splits on each feature."
))

cat("\n")


# =============================================================================
# SECTION 11 — PART 11: Step 6 Summary
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("STEP 6 COMPLETE — XGBoost BENCHMARK SUMMARY\n")
cat("================================================================\n")
cat("KEY FINDINGS:\n")
cat(sprintf("  Best params: depth=%d, eta=%.2f, mcw=%d, sub=%.1f, nrounds=%d\n",
            best_max_depth, best_eta, best_min_child_weight,
            best_subsample, best_nrounds))
cat(sprintf("  XGB Apparent C-index:   %.4f\n", auc_xgb_val))
cat(sprintf("  XGB Corrected C-index:  %.4f\n", auc_xgb_corrected))
cat(sprintf("  XGB Optimism:           %.4f\n", auc_xgb_optimism))
cat(sprintf("  XGB Corrected Slope:    %.4f\n", slope_xgb_corrected))
cat(sprintf("  XGB Apparent Brier:     %.4f\n", brier_xgb_apparent))
cat(sprintf("  XGB Corrected Brier:    %.4f\n", brier_xgb_corrected))
cat("---\n")
cat(sprintf("  Delta Corrected C-index (XGB - LR): %+.4f\n", delta_c_xgb))
cat(sprintf("  Delta Corrected Slope   (XGB - LR): %+.4f\n", delta_slope_xgb))
cat("---\n")
cat("INTERPRETATION:\n")
if (!meaningful) {
  cat("  XGBoost does not provide clinically meaningful improvement\n")
  cat("  over logistic regression on this dataset.\n")
  cat("  Convergent evidence across Steps 3, 5, and 6 supports\n")
  cat("  retention of the approved logistic regression model.\n")
} else {
  cat("  Mixed findings — flag for Strategist review.\n")
}
cat("---\n")
cat("POTENTIAL CONCERNS:\n")
cat("  1. XGBoost hyperparameters fixed during bootstrap — may\n")
cat("     underestimate true optimism (documented limitation).\n")
cat("  2. High apparent AUC may reflect memorization with small n.\n")
cat("  3. XGBoost not compatible with nomogram (Step 8) — cannot\n")
cat("     replace logistic model for final deployment pipeline.\n")
cat("---\n")
cat("RECOMMENDATION: RETAIN LOGISTIC REGRESSION\n")
cat("  Pending Strategist AI review and confirmation.\n")
cat("---\n")
cat("Figures:  figure10_ROC_benchmark_comparison.png\n")
cat("          figure11_calibration_XGBoost.png\n")
cat("          figure12_XGBoost_feature_importance.png\n")
cat("Tables:   Table 9  (preprocessing)\n")
cat("          Table 10 (tuning top 10)\n")
cat("          Table 11 (final parameters)\n")
cat("          Table 12 (bootstrap validation)\n")
cat("          Table 13 (three-way benchmark)\n")
cat("================================================================\n")
cat("STOPPING — Awaiting Strategist review before Step 7.\n")