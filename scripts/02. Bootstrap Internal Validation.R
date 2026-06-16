# =============================================================================
# PROJECT: Clinical Prediction Model — Occult Distant Metastasis in Gastric Cancer
# STEP 2:  Bootstrap Internal Validation
# AUTHOR:  Nhan Tran
# DATE:    160626
# TRIPOD PHASE: Internal validation via bootstrap (B = 2000)
# DEPENDENCY: Step 1 must be run first (model_glm, model_lrm, df, dd must exist)
# =============================================================================


# =============================================================================
# SECTION 0: Environment & Reproducibility
# =============================================================================

set.seed(2025)  # Same seed as Step 1 for reproducibility chain

# Required packages (should already be loaded from Step 1)
library(rms)       # validate(), calibrate(), lrm()
library(ggplot2)   # Publication figures
library(dplyr)     # Data wrangling
library(gt)        # Publication tables

# Bootstrap iterations
B_validate  <- 2000   # For validate() — optimism correction
B_calibrate <- 2000   # For calibrate() — calibration curve

# If your machine is slow, reduce both to 1000 and note this in your manuscript


# =============================================================================
# SECTION 1: Confirm Step 1 Objects Are Available
# =============================================================================

# These objects MUST exist from Step 1:
#   df           — full dataset (n = 114)
#   model_glm    — glm() fitted model
#   model_lrm    — lrm() fitted model (x = TRUE, y = TRUE)
#   dd           — datadist object
#   pred_probs   — predicted probabilities from glm (for Brier)

# Safety checks
if (!exists("model_lrm")) stop("model_lrm not found. Run Step 1 first.")
if (!exists("model_glm")) stop("model_glm not found. Run Step 1 first.")
if (!exists("df"))        stop("df not found. Run Step 1 first.")

cat("=== STEP 2 PREREQUISITES CHECK ===\n")
cat("model_lrm:  OK\n")
cat("model_glm:  OK\n")
cat("df rows:    ", nrow(df), "\n")
cat("Events (M1):", sum(df$Metastasis_TNM), "\n")
cat("Non-events: ", sum(df$Metastasis_TNM == 0), "\n")
cat("Event rate: ", round(mean(df$Metastasis_TNM), 3), "\n\n")

# Confirm datadist is still set
options(datadist = "dd")


# =============================================================================
# SECTION 2 — PART 1: Optimism-Corrected Validation via validate()
#
# rms::validate() implements the Harrell bootstrap optimism correction.
# It refits the model B times on bootstrap samples and evaluates on original data.
# Metrics returned: Dxy, R2, Intercept, Slope, Emax, D, U, Q, Brier
# =============================================================================

cat("=== PART 1: BOOTSTRAP OPTIMISM CORRECTION ===\n")
cat("Running validate() with B =", B_validate, "bootstrap iterations...\n")
cat("(This may take 30–90 seconds)\n\n")

val_results <- validate(
  model_lrm,
  method = "boot",   # Bootstrap resampling (Harrell method)
  B      = B_validate
)

# Print raw validate() output
print(val_results)

# --- Extract key metrics ---

# Dxy (Somers' D) — apparent, optimism, corrected
dxy_apparent  <- val_results["Dxy", "index.orig"]
dxy_optimism  <- val_results["Dxy", "optimism"]
dxy_corrected <- val_results["Dxy", "index.corrected"]

# C-index = Dxy/2 + 0.5
cindex_apparent  <- dxy_apparent  / 2 + 0.5
cindex_optimism  <- dxy_optimism  / 2          # Optimism on C-index scale
cindex_corrected <- dxy_corrected / 2 + 0.5

# Calibration slope — apparent, optimism, corrected
slope_apparent  <- val_results["Slope", "index.orig"]
slope_optimism  <- val_results["Slope", "optimism"]
slope_corrected <- val_results["Slope", "index.corrected"]

# Calibration intercept
intercept_apparent  <- val_results["Intercept", "index.orig"]
intercept_optimism  <- val_results["Intercept", "optimism"]
intercept_corrected <- val_results["Intercept", "index.corrected"]

# Brier Score from validate() (note: rms scales differently — see Part 5)
brier_val_apparent  <- val_results["B",  "index.orig"]
brier_val_corrected <- val_results["B",  "index.corrected"]

cat("\n--- EXTRACTED METRICS ---\n")
cat(sprintf("C-index  | Apparent: %.4f | Optimism: %.4f | Corrected: %.4f\n",
            cindex_apparent, cindex_optimism, cindex_corrected))
cat(sprintf("Slope    | Apparent: %.4f | Optimism: %.4f | Corrected: %.4f\n",
            slope_apparent, slope_optimism, slope_corrected))
cat(sprintf("Intercept| Apparent: %.4f | Optimism: %.4f | Corrected: %.4f\n",
            intercept_apparent, intercept_optimism, intercept_corrected))


# =============================================================================
# SECTION 3 — PART 2: Apparent vs Corrected Comparison Table
#
# NOTE: For logistic regression, AUC = C-index
# This is stated explicitly in the table footnote per TRIPOD convention
# =============================================================================

cat("\n=== PART 2: COMPARISON TABLE ===\n")

# Apparent Brier from Step 1 (glm-based, standard formula)
if (!exists("pred_probs")) {
  pred_probs <- predict(model_glm, type = "response")
}
brier_apparent <- mean((pred_probs - df$Metastasis_TNM)^2)

# Build comparison data frame
comparison_table <- data.frame(
  Metric = c(
    "C-index (= AUC)",
    "Somers' Dxy",
    "Calibration Slope",
    "Calibration Intercept",
    "Brier Score"
  ),
  Apparent = c(
    round(cindex_apparent,    3),
    round(dxy_apparent,       3),
    round(slope_apparent,     3),
    round(intercept_apparent, 3),
    round(brier_apparent,     4)
  ),
  Optimism = c(
    round(cindex_optimism,    3),
    round(dxy_optimism,       3),
    round(slope_optimism,     3),
    round(intercept_optimism, 3),
    NA    # Brier optimism from manual loop in Part 5
  ),
  Corrected = c(
    round(cindex_corrected,    3),
    round(dxy_corrected,       3),
    round(slope_corrected,     3),
    round(intercept_corrected, 3),
    NA    # To be filled after Part 5
  )
)

print(comparison_table)

# gt publication table
comparison_gt <- comparison_table %>%
  gt() %>%
  tab_header(
    title    = "Table 2. Bootstrap Internal Validation Results",
    subtitle = paste0("Baseline logistic regression model (B = ", B_validate,
                      " bootstrap iterations)")
  ) %>%
  cols_label(
    Metric    = "Performance Metric",
    Apparent  = "Apparent",
    Optimism  = "Optimism",
    Corrected = "Optimism-Corrected"
  ) %>%
  sub_missing(missing_text = "See Part 5") %>%
  tab_footnote(
    footnote = paste0(
      "C-index = AUC for binary logistic regression. ",
      "Optimism estimated via ", B_validate, "-iteration bootstrap (Harrell method). ",
      "Corrected = Apparent − Optimism. ",
      "Calibration slope: ideal = 1.0. ",
      "Calibration intercept: ideal = 0.0."
    )
  ) %>%
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_column_labels()
  )

print(comparison_gt)

# Optional: save
# gtsave(comparison_gt, "table2_bootstrap_validation.html")


# =============================================================================
# SECTION 4 — PART 3: Bootstrap Calibration Curve
#
# calibrate() resamples and compares predicted vs observed probabilities
# Produces: apparent curve, bias-corrected curve, ideal diagonal
# =============================================================================

cat("\n=== PART 3: BOOTSTRAP CALIBRATION CURVE ===\n")
cat("Running calibrate() with B =", B_calibrate, "iterations...\n")

cal_boot <- calibrate(
  model_lrm,
  method = "boot",
  B      = B_calibrate
)

# --- Extract calibration data for ggplot2 ---
# Use actual column names from calibrate() output
print(colnames(cal_boot))  # Confirm names before extracting

cal_df <- as.data.frame(cal_boot[, c("predy",
                                     "calibrated.orig",
                                     "calibrated.corrected")])
colnames(cal_df) <- c("Predicted", "Apparent", "Corrected")
cat("\nCalibration data (first 6 rows):\n")
print(head(cal_df))

# --- Publication-ready calibration figure ---
cal_plot <- ggplot(cal_df, aes(x = Predicted)) +
  
  # Ideal line
  geom_abline(slope     = 1,
              intercept = 0,
              linetype  = "solid",
              color     = "gray60",
              linewidth = 0.8) +
  
  # Apparent calibration
  geom_line(aes(y = Apparent, color = "Apparent"),
            linewidth = 1.0, linetype = "dashed") +
  geom_point(aes(y = Apparent, color = "Apparent"),
             size = 2.5) +
  
  # Bias-corrected calibration
  geom_line(aes(y = Corrected, color = "Bias-corrected"),
            linewidth = 1.2) +
  geom_point(aes(y = Corrected, color = "Bias-corrected"),
             size = 2.5) +
  
  # Scales and labels
  scale_color_manual(
    name   = "Calibration",
    values = c("Apparent"       = "#E07B39",
               "Bias-corrected" = "#2C3E8C")
  ) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  
  labs(
    title    = "Calibration Plot — Baseline Logistic Regression",
    subtitle = paste0("Bootstrap internal validation (B = ", B_calibrate, ")"),
    x        = "Predicted Probability",
    y        = "Observed Proportion (Fraction of Events)",
    caption  = "Gray diagonal = ideal calibration"
  ) +
  
  theme_classic(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold"),
    plot.subtitle = element_text(color = "gray40"),
    legend.position = c(0.15, 0.85),
    legend.background = element_rect(fill = "white", color = "gray80")
  )

print(cal_plot)

ggsave("figure3_calibration_bootstrap.png", cal_plot,
       width = 6, height = 5.5, dpi = 300)

cat("Calibration figure saved: figure3_calibration_bootstrap.png\n")


# =============================================================================
# SECTION 5 — PART 4: Written Overfitting Assessment
#
# Thresholds used (standard prediction modeling conventions):
#   Optimism in C-index < 0.01  → negligible
#                       0.01–0.03 → small
#                       > 0.03    → moderate/large
#   Calibration slope corrected > 0.90 → acceptable
#                               0.80–0.90 → moderate shrinkage needed
#                               < 0.80   → substantial overfitting
# =============================================================================

cat("\n=== PART 4: OVERFITTING ASSESSMENT ===\n\n")

# --- Overfitting severity ---
optimism_level <- ifelse(
  cindex_optimism < 0.01, "NEGLIGIBLE",
  ifelse(cindex_optimism < 0.03, "SMALL", "MODERATE TO LARGE")
)

# --- Calibration slope interpretation ---
slope_interp <- ifelse(
  slope_corrected > 0.90, "ACCEPTABLE (>0.90)",
  ifelse(slope_corrected > 0.80,
         "MODERATE SHRINKAGE INDICATED (0.80–0.90)",
         "SUBSTANTIAL OVERFITTING (<0.80)")
)

cat("--- DISCRIMINATION ---\n")
cat(sprintf("  Apparent C-index:   %.3f\n", cindex_apparent))
cat(sprintf("  Optimism:           %.3f  [%s]\n", cindex_optimism, optimism_level))
cat(sprintf("  Corrected C-index:  %.3f\n", cindex_corrected))
cat("\n")

cat("--- CALIBRATION ---\n")
cat(sprintf("  Corrected Slope:     %.3f  [%s]\n", slope_corrected, slope_interp))
cat(sprintf("  Corrected Intercept: %.3f  (ideal = 0.00)\n", intercept_corrected))
cat("\n")

cat("--- GENERALIZABILITY ASSESSMENT ---\n")
if (cindex_optimism < 0.03 && slope_corrected > 0.85) {
  cat("  CONCLUSION: Model shows LIMITED evidence of overfitting.\n")
  cat("  The corrected C-index suggests acceptable discriminative ability\n")
  cat("  in future patients drawn from the same population.\n")
  cat("  Calibration slope near 1.0 indicates predictions are not\n")
  cat("  excessively extreme. Model likely generalizes reasonably\n")
  cat("  within this clinical context.\n")
} else if (cindex_optimism >= 0.03 && slope_corrected > 0.80) {
  cat("  CONCLUSION: Model shows MODERATE overfitting.\n")
  cat("  Optimism in C-index exceeds 0.03, suggesting the apparent\n")
  cat("  performance is meaningfully inflated. Calibration slope <1.0\n")
  cat("  indicates predictions may be too extreme.\n")
  cat("  Consider: penalized regression (Elastic Net, Step 5) as benchmark.\n")
  cat("  Shrinkage of predictions may improve calibration on new data.\n")
} else {
  cat("  CONCLUSION: Model shows evidence of SUBSTANTIAL overfitting.\n")
  cat("  Both discrimination optimism and calibration slope suggest\n")
  cat("  considerable performance inflation in apparent metrics.\n")
  cat("  Penalized or shrinkage methods are strongly recommended.\n")
}

cat("\n  NOTE: With n=114 and 51 events (EPV ≈ 17 per predictor term),\n")
cat("  this sample size is moderate for the number of predictors.\n")
cat("  External validation remains necessary before clinical deployment.\n")


# =============================================================================
# SECTION 6 — PART 5: Bootstrap-Corrected Brier Score
#
# rms::validate() reports a "B" (Brier-like) row but uses a scaled version.
# We compute the standard Brier Score manually via bootstrap loop
# to remain consistent with Step 1's Brier calculation.
#
# METHOD:
#   For each bootstrap b:
#     1. Fit model on bootstrap sample
#     2. Get predicted probs on bootstrap sample → Brier_boot
#     3. Get predicted probs on original data    → Brier_orig
#     4. optimism_b = Brier_boot - Brier_orig
#   Corrected Brier = Apparent Brier - mean(optimism_b)
# =============================================================================

cat("\n=== PART 5: BOOTSTRAP-CORRECTED BRIER SCORE ===\n")
cat("Running manual bootstrap loop for Brier Score correction...\n")
cat("B = ", B_validate, "iterations\n\n")

set.seed(2025)
n <- nrow(df)

brier_optimism_vec <- numeric(B_validate)

for (b in seq_len(B_validate)) {
  
  # Step 1: Bootstrap sample indices
  boot_idx  <- sample(seq_len(n), size = n, replace = TRUE)
  boot_data <- df[boot_idx, ]
  
  # Step 2: Refit glm on bootstrap sample
  # We use glm here (not lrm) for speed; both will give identical Brier estimates
  boot_model <- tryCatch(
    glm(
      Metastasis_TNM ~ Tumor_mlocation3 + NLR + PLR,
      data   = boot_data,
      family = binomial(link = "logit")
    ),
    error = function(e) NULL   # Skip failed fits (rare, e.g., separation)
  )
  
  if (is.null(boot_model)) next
  
  # Step 3: Brier on bootstrap sample (apparent performance in bootstrap world)
  pred_boot <- predict(boot_model, newdata = boot_data, type = "response")
  brier_boot <- mean((pred_boot - boot_data$Metastasis_TNM)^2)
  
  # Step 4: Brier on original data (test performance)
  pred_orig <- predict(boot_model, newdata = df, type = "response")
  brier_orig <- mean((pred_orig - df$Metastasis_TNM)^2)
  
  # Step 5: Optimism for this iteration
  brier_optimism_vec[b] <- brier_boot - brier_orig
}

# Mean optimism across all iterations
brier_optimism_mean <- mean(brier_optimism_vec, na.rm = TRUE)

# Corrected Brier = Apparent - Optimism
brier_corrected <- brier_apparent - brier_optimism_mean

# Null Brier (reference: always predict prevalence)
prevalence  <- mean(df$Metastasis_TNM)
brier_null  <- mean((prevalence - df$Metastasis_TNM)^2)

# Scaled Brier (proportion of null explained, higher = better)
scaled_brier_apparent  <- 1 - (brier_apparent  / brier_null)
scaled_brier_corrected <- 1 - (brier_corrected / brier_null)

cat(sprintf("Apparent Brier Score:          %.4f\n", brier_apparent))
cat(sprintf("Optimism (bootstrap):          %.4f\n", brier_optimism_mean))
cat(sprintf("Corrected Brier Score:         %.4f\n", brier_corrected))
cat(sprintf("Null Brier Score:              %.4f\n", brier_null))
cat(sprintf("Scaled Brier (Apparent):       %.4f\n", scaled_brier_apparent))
cat(sprintf("Scaled Brier (Corrected):      %.4f\n", scaled_brier_corrected))

# --- Update comparison table with Brier values ---
comparison_table$Optimism[5]  <- round(brier_optimism_mean, 4)
comparison_table$Corrected[5] <- round(brier_corrected, 4)

cat("\n=== UPDATED COMPARISON TABLE (with Brier) ===\n")
print(comparison_table)


# =============================================================================
# SECTION 7: Step 2 Summary
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("STEP 2 COMPLETE — BOOTSTRAP INTERNAL VALIDATION SUMMARY\n")
cat("================================================================\n")
cat(sprintf("Bootstrap iterations:          B = %d\n", B_validate))
cat(sprintf("Dataset:                       n = %d, Events = %d\n",
            nrow(df), sum(df$Metastasis_TNM)))
cat("---\n")
cat(sprintf("C-index    | Apparent: %.3f | Optimism: %.3f | Corrected: %.3f\n",
            cindex_apparent, cindex_optimism, cindex_corrected))
cat(sprintf("Cal. Slope | Apparent: %.3f | Optimism: %.3f | Corrected: %.3f\n",
            slope_apparent, slope_optimism, slope_corrected))
cat(sprintf("Brier      | Apparent: %.4f | Optimism: %.4f | Corrected: %.4f\n",
            brier_apparent, brier_optimism_mean, brier_corrected))
cat("---\n")
cat("Figures:    figure3_calibration_bootstrap.png\n")
cat("Tables:     Table 2 (comparison_gt)\n")
cat("================================================================\n")
cat("STOPPING — Awaiting review before Step 3 (Restricted Cubic Splines).\n")