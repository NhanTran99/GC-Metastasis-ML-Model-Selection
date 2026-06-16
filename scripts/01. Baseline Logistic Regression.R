# =============================================================================
# PROJECT: Clinical Prediction Model — Occult Distant Metastasis in Gastric Cancer
# STEP 1:  Baseline Logistic Regression (Apparent Performance)
# AUTHOR:  Nhan Tran
# DATE:    150626
# TRIPOD PHASE: Development (full dataset, n = 114)
# =============================================================================

# =============================================================================
# SECTION 0: Environment Setup
# =============================================================================

# Set seed for reproducibility (relevant for any jitter in plots)
set.seed(2025)

# Required packages
# Install if not already present:
# install.packages(c("rms", "pROC", "ggplot2", "gt", "dplyr", "CalibrationCurves"))

library(rms)       # lrm(), calibrate(), nomogram()
library(pROC)      # roc(), auc(), ci.auc()
library(ggplot2)   # publication-ready figures
library(dplyr)     # data wrangling
library(gt)        # publication-ready tables


# =============================================================================
# SECTION 1: Data Loading & Integrity Checks
# =============================================================================

# --- Load your dataset ---
# Replace this line with your actual file path:
df <- read.delim("NLR.txt", sep = "\t", header = TRUE)

# For now, a placeholder check structure is defined.
# Uncomment and replace with your actual import method.

# df <- read.csv("gastric_metastasis_data.csv")

# --- Integrity checks ---
# Run these after loading to confirm dataset matches specification

cat("=== DATASET INTEGRITY CHECK ===\n")
cat("N rows:", nrow(df), "\n")           # Expected: 114
cat("N cols:", ncol(df), "\n")
cat("Missing values:\n")
print(colSums(is.na(df)))               # Expected: all 0

# Confirm outcome variable
cat("\nOutcome distribution (Metastasis_TNM):\n")
print(table(df$Metastasis_TNM))         # Should be 0/1 binary

# Confirm predictor types
cat("\nPredictor classes:\n")
cat("Tumor_mlocation3:", class(df$Tumor_mlocation3), "\n")
cat("NLR:", class(df$NLR), "\n")
cat("PLR:", class(df$PLR), "\n")

# IMPORTANT: Tumor_mlocation3 must be a factor for correct dummy coding
# If it is stored as numeric/character, convert here:
df$Tumor_mlocation3 <- factor(df$Tumor_mlocation3)

# Confirm factor levels (reference category will be the first level)
cat("\nTumor_mlocation3 levels (first = reference):\n")
print(levels(df$Tumor_mlocation3))


# =============================================================================
# SECTION 2: Model Fitting
# =============================================================================

# --- 2A: glm() fit ---
# Used for: coefficient table, prediction function, Brier score
# Standard binomial logistic regression, no penalization

model_glm <- glm(
  Metastasis_TNM ~ Tumor_mlocation3 + NLR + PLR,
  data   = df,
  family = binomial(link = "logit")
)

cat("\n=== GLM MODEL SUMMARY ===\n")
print(summary(model_glm))


# --- 2B: rms::lrm() fit ---
# Used for: calibration plot (calibrate()), future RCS (Step 3), nomogram (Step 8)
# Requires rms datadist for downstream functions

dd <- datadist(df)
options(datadist = "dd")

model_lrm <- lrm(
  Metastasis_TNM ~ Tumor_mlocation3 + NLR + PLR,
  data  = df,
  x     = TRUE,   # Store design matrix (required for calibrate())
  y     = TRUE    # Store outcome vector (required for calibrate())
)

cat("\n=== LRM MODEL SUMMARY ===\n")
print(model_lrm)


# =============================================================================
# SECTION 3: Publication-Ready Coefficient Table
# Reports: Coefficient, SE, OR, 95% CI, p-value
# =============================================================================

# Extract from glm object
coef_table <- data.frame(
  Term       = names(coef(model_glm)),
  Coefficient = round(coef(model_glm), 4),
  SE         = round(summary(model_glm)$coefficients[, "Std. Error"], 4),
  OR         = round(exp(coef(model_glm)), 3),
  CI_lower   = round(exp(confint(model_glm))[, 1], 3),
  CI_upper   = round(exp(confint(model_glm))[, 2], 3),
  p_value    = round(summary(model_glm)$coefficients[, "Pr(>|z|)"], 4)
)

# Format OR (95% CI) as single string for publication
coef_table <- coef_table %>%
  mutate(
    `OR (95% CI)` = paste0(OR, " (", CI_lower, "–", CI_upper, ")")
  ) %>%
  select(Term, Coefficient, SE, `OR (95% CI)`, p_value)

cat("\n=== PUBLICATION COEFFICIENT TABLE ===\n")
print(coef_table)

# Render as gt table (export-ready)
coef_gt <- coef_table %>%
  gt() %>%
  tab_header(
    title    = "Table 1. Logistic Regression Model Coefficients",
    subtitle = "Outcome: Occult Distant Metastasis (Metastasis_TNM)"
  ) %>%
  cols_label(
    Term        = "Variable",
    Coefficient = "β Coefficient",
    SE          = "Std. Error",
    `OR (95% CI)` = "OR (95% CI)",
    p_value     = "p-value"
  ) %>%
  fmt_number(columns = c(Coefficient, SE), decimals = 3) %>%
  fmt_number(columns = p_value, decimals = 4) %>%
  tab_footnote(
    footnote = "OR = Odds Ratio; CI = Confidence Interval; Reference category for Tumor_mlocation3 = first factor level."
  )

print(coef_gt)

# Optional: save as HTML
# gtsave(coef_gt, "table1_coefficients.html")


# =============================================================================
# SECTION 4: ROC Curve & AUC with 95% CI
# =============================================================================

# Predicted probabilities from glm
pred_probs <- predict(model_glm, type = "response")

# ROC object
roc_obj <- roc(
  response  = df$Metastasis_TNM,
  predictor = pred_probs,
  ci        = TRUE,
  direction = "<"     # Higher probability = higher risk of metastasis
)

# AUC and 95% CI
auc_val    <- round(as.numeric(auc(roc_obj)), 3)
auc_ci     <- round(as.numeric(ci.auc(roc_obj)), 3)

cat("\n=== DISCRIMINATION ===\n")
cat("AUC:", auc_val, "\n")
cat("95% CI:", auc_ci[1], "–", auc_ci[3], "\n")

# --- ROC Plot (ggplot2, publication-ready) ---
roc_df <- data.frame(
  Specificity = roc_obj$specificities,
  Sensitivity = roc_obj$sensitivities
)

roc_plot <- ggplot(roc_df, aes(x = 1 - Specificity, y = Sensitivity)) +
  geom_line(color = "#2C3E8C", linewidth = 1.2) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "gray50") +
  annotate("text",
           x    = 0.65, y = 0.15,
           label = paste0("AUC = ", auc_val,
                          "\n95% CI: ", auc_ci[1], "–", auc_ci[3]),
           size = 4, hjust = 0) +
  labs(
    title    = "ROC Curve — Baseline Logistic Regression",
    subtitle = "Apparent performance (full dataset, n = 114)",
    x        = "1 − Specificity (False Positive Rate)",
    y        = "Sensitivity (True Positive Rate)"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold"),
    plot.subtitle = element_text(color = "gray40")
  )

print(roc_plot)

# Save figure
ggsave("figure1_ROC_baseline.png", roc_plot,
       width = 6, height = 5, dpi = 300)

cat("ROC figure saved: figure1_ROC_baseline.png\n")


# =============================================================================
# SECTION 5: Apparent Calibration Plot
# Uses rms::calibrate() — loess-smoothed calibration on apparent data
# NOTE: This is apparent calibration only. Bootstrap correction in Step 2.
# =============================================================================

# calibrate() with B = 0 gives apparent calibration (no resampling)
cal_apparent <- calibrate(
  model_lrm,
  method = "boot",
  B      = 1
)

png("figure2_calibration_apparent.png", width = 700, height = 600, res = 150)
plot(cal_apparent,
     main = "Calibration Plot — Baseline Model (Apparent)",
     xlab = "Predicted Probability",
     ylab = "Observed Proportion")
dev.off()


# Base R calibration plot (rms default)
cat("\n=== CALIBRATION PLOT ===\n")
cat("Plotting apparent calibration (rms::calibrate)...\n")

png("figure2_calibration_apparent.png", width = 700, height = 600, res = 150)
plot(cal_apparent,
     main = "Calibration Plot — Baseline Model (Apparent)",
     xlab = "Predicted Probability",
     ylab = "Observed Proportion")
dev.off()

cat("Calibration figure saved: figure2_calibration_apparent.png\n")


# =============================================================================
# SECTION 6: Brier Score (Overall Performance)
# Brier Score = mean squared error of predicted probabilities
# Range: 0 (perfect) to 1 (worst); uninformative model ≈ prevalence*(1-prevalence)
# =============================================================================

brier_score <- mean((pred_probs - df$Metastasis_TNM)^2)

# Reference Brier score (null model: predict prevalence for everyone)
prevalence        <- mean(df$Metastasis_TNM)
brier_null        <- mean((prevalence - df$Metastasis_TNM)^2)
brier_scaled      <- 1 - (brier_score / brier_null)   # Scaled Brier (higher = better)

cat("\n=== OVERALL PERFORMANCE ===\n")
cat("Brier Score (apparent):", round(brier_score, 4), "\n")
cat("Null Brier Score:       ", round(brier_null, 4), "\n")
cat("Scaled Brier Score:     ", round(brier_scaled, 4),
    "(proportion of null explained)\n")


# =============================================================================
# SECTION 7: Reusable Prediction Function
# Accepts new patient values and returns predicted probability
# Designed for reuse in Step 10 (web calculator)
# =============================================================================

predict_metastasis <- function(Tumor_mlocation3,
                               NLR,
                               PLR,
                               model = model_glm) {
  # Input validation
  stopifnot(is.numeric(NLR), is.numeric(PLR))
  stopifnot(length(NLR) == length(PLR))
  
  # Build new data frame
  new_data <- data.frame(
    Tumor_mlocation3 = factor(Tumor_mlocation3,
                              levels = levels(df$Tumor_mlocation3)),
    NLR              = NLR,
    PLR              = PLR
  )
  
  # Return predicted probability
  prob <- predict(model, newdata = new_data, type = "response")
  return(round(prob, 4))
}

# --- Example usage ---
cat("\n=== PREDICTION FUNCTION TEST ===\n")
test_prob <- predict_metastasis(
  Tumor_mlocation3 = levels(df$Tumor_mlocation3)[1],  # Reference level
  NLR              = 3.5,
  PLR              = 150
)
cat("Example predicted probability:", test_prob, "\n")


# =============================================================================
# SECTION 8: Step 1 Summary Output
# =============================================================================

cat("\n")
cat("========================================================\n")
cat("STEP 1 COMPLETE — BASELINE MODEL SUMMARY\n")
cat("========================================================\n")
cat("Model:        glm + lrm (Metastasis_TNM ~ Tumor_mlocation3 + NLR + PLR)\n")
cat("N:            114\n")
cat("AUC:          ", auc_val, "(95% CI:", auc_ci[1], "–", auc_ci[3], ")\n")
cat("Brier Score:  ", round(brier_score, 4), "\n")
cat("Scaled Brier: ", round(brier_scaled, 4), "\n")
cat("Figures:      figure1_ROC_baseline.png\n")
cat("              figure2_calibration_apparent.png\n")
cat("========================================================\n")
cat("STOPPING — Awaiting review before Step 2.\n")
