# =============================================================================
# PROJECT: Clinical Prediction Model — Occult Distant Metastasis in Gastric Cancer
# STEP 4:  Decision Curve Analysis (DCA)
# AUTHOR:  Nhan Tran
# DATE:    160626
# DEPENDENCY: Steps 1-3 must be run first.
#             Required objects: df, model_glm
# =============================================================================


# =============================================================================
# SECTION 0: Environment & Prerequisites
# =============================================================================

set.seed(2025)

library(ggplot2)
library(dplyr)
library(gt)

# --- Install/load DCA package ---
# Primary: dcurves (modern, ggplot2-native, CRAN-stable)
# Fallback: rmda (older, base-R plots)
#
# We use dcurves as primary — it is actively maintained, handles binary
# outcomes directly, and produces ggplot2 output for publication formatting.
#
# install.packages("dcurves")
library(dcurves)

# Safety checks
if (!exists("model_glm")) stop("model_glm not found. Run Steps 1-3 first.")
if (!exists("df"))        stop("df not found. Run Steps 1-3 first.")

cat("=== STEP 4 PREREQUISITES CHECK ===\n")
cat("n =", nrow(df), "\n")
cat("Events (M1) =", sum(df$Metastasis_TNM), "\n")
cat("Prevalence =", round(mean(df$Metastasis_TNM), 3), "\n\n")


# =============================================================================
# SECTION 1 — PART 1: Predicted Probabilities & Summary Statistics
# =============================================================================

cat("=== PART 1: PREDICTED PROBABILITIES ===\n")

# Generate predicted probabilities from approved linear model
df$pred_prob <- predict(model_glm, type = "response")

# Summary statistics
prob_summary <- data.frame(
  Statistic = c("Minimum", "25th Percentile", "Median",
                "Mean", "75th Percentile", "Maximum"),
  Value = round(c(
    min(df$pred_prob),
    quantile(df$pred_prob, 0.25),
    median(df$pred_prob),
    mean(df$pred_prob),
    quantile(df$pred_prob, 0.75),
    max(df$pred_prob)
  ), 4)
)

print(prob_summary)

# Distribution by outcome group — clinical sanity check
cat("\nMean predicted probability by outcome:\n")
df %>%
  group_by(Metastasis_TNM) %>%
  summarise(
    n         = n(),
    mean_prob = round(mean(pred_prob), 3),
    sd_prob   = round(sd(pred_prob), 3),
    .groups   = "drop"
  ) %>%
  print()

# Histogram of predicted probabilities
p_hist <- ggplot(df, aes(x = pred_prob,
                         fill = factor(Metastasis_TNM))) +
  geom_histogram(binwidth = 0.05, position = "identity",
                 alpha = 0.6, color = "white") +
  scale_fill_manual(
    values = c("0" = "#4A90D9", "1" = "#E07B39"),
    labels = c("0" = "M0 (no metastasis)", "1" = "M1 (metastasis)"),
    name   = "Outcome"
  ) +
  labs(
    title    = "Distribution of Predicted Probabilities",
    subtitle = "Approved linear logistic regression model",
    x        = "Predicted Probability of Distant Metastasis",
    y        = "Count"
  ) +
  theme_classic(base_size = 13) +
  theme(plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(color = "gray40"))

print(p_hist)
ggsave("figure6_predicted_prob_distribution.png", p_hist,
       width = 6.5, height = 4.5, dpi = 300)

cat("Probability distribution figure saved.\n\n")


# =============================================================================
# SECTION 2 — PARTS 2 & 3: Decision Curve Analysis
#
# dcurves::dca() computes net benefit across a threshold range
# thresholds: seq(0.05, 0.50, by = 0.01) — 5% to 50%
# Strategies: treat_none (automatic), treat_all (automatic), nomogram
# =============================================================================

cat("=== PARTS 2 & 3: DECISION CURVE ANALYSIS ===\n")

# Run DCA
# dcurves uses the outcome and predictor columns directly from the dataframe
dca_result <- dca(
  formula    = Metastasis_TNM ~ pred_prob,
  data       = df,
  thresholds = seq(0.05, 0.50, by = 0.01)
)

# Print net benefit summary
cat("\nDCA result summary (first 10 rows):\n")
print(head(as.data.frame(dca_result$dca), 10))

# Fix: extract directly from dca_result without filtering
dca_df_full <- as.data.frame(dca_result$dca)

# Check what's actually in the data
cat("Columns in dca_result$dca:\n")
print(colnames(dca_df_full))
cat("\nUnique variable values:\n")
print(unique(dca_df_full$variable))


# =============================================================================
# SECTION 3 — PART 4: Publication-Quality DCA Figure
#
# Extract net benefit data from dca_result for manual ggplot2 construction
# This gives full control over appearance for manuscript submission
# =============================================================================

cat("\n=== PART 4: DCA FIGURE ===\n")

# Extract net benefit data
# Extract and recode strategy labels using correct variable names
dca_df <- dca_df_full %>%
  filter(variable %in% c("all", "none", "pred_prob")) %>%
  mutate(
    Strategy = case_when(
      variable == "all"       ~ "Treat All",
      variable == "none"      ~ "Treat None",
      variable == "pred_prob" ~ "Nomogram"
    )
  )


# Publication-ready DCA plot
dca_plot <- ggplot(dca_df,
                   aes(x     = threshold,
                       y     = net_benefit,
                       color = Strategy,
                       linetype = Strategy)) +
  
  geom_line(linewidth = 1.1, na.rm = TRUE) +
  
  # Reference line at net benefit = 0
  geom_hline(yintercept = 0,
             linetype   = "dotted",
             color      = "gray60",
             linewidth  = 0.6) +
  
  scale_color_manual(
    values = c(
      "Treat All"  = "#999999",
      "Treat None" = "#000000",
      "Nomogram"   = "#2C3E8C"
    )
  ) +
  scale_linetype_manual(
    values = c(
      "Treat All"  = "dashed",
      "Treat None" = "solid",
      "Nomogram"   = "solid"
    )
  ) +
  
  scale_x_continuous(
    limits = c(0.05, 0.50),
    breaks = seq(0.05, 0.50, by = 0.05),
    labels = scales::percent_format(accuracy = 1)
  ) +
  scale_y_continuous(
    breaks = seq(-0.1, 0.5, by = 0.05)
  ) +
  
  coord_cartesian(ylim = c(-0.05, 0.55)) +
  
  labs(
    title    = "Decision Curve Analysis",
    subtitle = paste0(
      "Nomogram for occult distant metastasis in gastric cancer\n",
      "Threshold range: 5%–50%"
    ),
    x        = "Threshold Probability",
    y        = "Net Benefit",
    color    = "Strategy",
    linetype = "Strategy",
    caption  = paste0(
      "Net benefit = (TP/n) − (FP/n) × (pt/(1−pt))\n",
      "n = 114; Events (M1) = 51; Prevalence = 44.7%"
    )
  ) +
  
  theme_classic(base_size = 13) +
  theme(
    plot.title      = element_text(face = "bold"),
    plot.subtitle   = element_text(color = "gray40", size = 10),
    plot.caption    = element_text(color = "gray50", size = 9),
    legend.position = c(0.80, 0.85),
    legend.background = element_rect(fill = "white", color = "gray80"),
    legend.title    = element_blank()
  )

print(dca_plot)
ggsave("figure7_DCA.png", dca_plot,
       width = 7, height = 5.5, dpi = 300)

cat("DCA figure saved: figure7_DCA.png\n")


# =============================================================================
# SECTION 4 — PART 5: Quantify Clinical Utility Threshold Region
#
# Identify thresholds where nomogram net benefit exceeds BOTH:
#   (a) Treat All
#   (b) Treat None (= 0)
# =============================================================================

cat("\n=== PART 5: THRESHOLD REGION OF CLINICAL UTILITY ===\n")

# Pivot to wide format for comparison
dca_wide <- dca_df %>%
  select(threshold, Strategy, net_benefit) %>%
  tidyr::pivot_wider(names_from  = Strategy,
                     values_from = net_benefit) %>%
  rename(
    nb_all      = `Treat All`,
    nb_none     = `Treat None`,
    nb_nomogram = Nomogram
  )

# Identify thresholds where nomogram beats both defaults
utility_region <- dca_wide %>%
  filter(
    !is.na(nb_nomogram),
    nb_nomogram > nb_all,
    nb_nomogram > nb_none,
    nb_nomogram > 0
  )

if (nrow(utility_region) > 0) {
  threshold_lower <- min(utility_region$threshold)
  threshold_upper <- max(utility_region$threshold)
  
  cat(sprintf(
    "Nomogram superior to both defaults: %.0f%% to %.0f%%\n",
    threshold_lower * 100,
    threshold_upper * 100
  ))
  
  # Net benefit at key thresholds
  cat("\nNet benefit at selected thresholds:\n")
  key_thresholds <- dca_wide %>%
    filter(threshold %in% c(0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40))
  
  print(
    key_thresholds %>%
      mutate(across(where(is.numeric), ~ round(.x, 4)))
  )
  
} else {
  threshold_lower <- NA
  threshold_upper <- NA
  cat("WARNING: Nomogram does not demonstrate superior net benefit\n")
  cat("over both default strategies in the 5%-50% range.\n")
  cat("Review DCA plot carefully before interpretation.\n")
}

# Net benefit summary table for publication
nb_table <- dca_wide %>%
  filter(threshold %in% seq(0.05, 0.50, by = 0.05)) %>%
  mutate(across(where(is.numeric), ~ round(.x, 4))) %>%
  select(threshold, nb_nomogram, nb_all, nb_none) %>%
  rename(
    `Threshold` = threshold,
    `Nomogram`  = nb_nomogram,
    `Treat All` = nb_all,
    `Treat None`= nb_none
  ) %>%
  mutate(`Threshold` = scales::percent(`Threshold`, accuracy = 1))

nb_gt <- nb_table %>%
  gt() %>%
  tab_header(
    title    = "Table 5. Net Benefit by Decision Threshold",
    subtitle = "Decision Curve Analysis — Approved Logistic Regression Model"
  ) %>%
  tab_footnote(
    footnote = paste0(
      "Net benefit computed at threshold probabilities from 5% to 50%. ",
      "Higher net benefit = greater clinical utility at that threshold. ",
      "Treat None net benefit = 0 by definition."
    )
  ) %>%
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_column_labels()
  )

print(nb_gt)


# =============================================================================
# SECTION 5 — PART 6: Clinical Interpretation
# =============================================================================

cat("\n=== PART 6: CLINICAL INTERPRETATION ===\n\n")

cat("CLINICAL CONTEXT:\n")
cat("  The model is intended to support prioritization of CT/MRI\n")
cat("  (not replace it) in patients evaluated by endoscopy + ultrasound + CBC.\n\n")

cat("INTERPRETATION FRAMEWORK:\n")
cat("  1. CLINICAL UTILITY:\n")
if (!is.na(threshold_lower)) {
  cat(sprintf(
    "     The nomogram demonstrates net benefit superior to both\n"
  ))
  cat(sprintf(
    "     treat-all and treat-none strategies across the threshold\n"
  ))
  cat(sprintf(
    "     range of approximately %.0f%%–%.0f%%.\n",
    threshold_lower * 100, threshold_upper * 100
  ))
  cat("     This suggests the model provides clinically actionable\n")
  cat("     risk stratification within this probability range.\n\n")
} else {
  cat("     The nomogram does not consistently outperform default\n")
  cat("     strategies. Clinical utility is limited.\n\n")
}

cat("  2. THRESHOLD INTERPRETATION:\n")
cat("     At lower thresholds (e.g., 10-15%): The nomogram may help\n")
cat("     identify which patients with seemingly low-risk profiles\n")
cat("     still warrant advanced imaging, reducing under-referral.\n\n")
cat("     At higher thresholds (e.g., 30-40%): The nomogram may\n")
cat("     help avoid unnecessary CT/MRI in patients whose predicted\n")
cat("     risk falls below the clinician's action threshold.\n\n")

cat("  3. UNNECESSARY IMAGING REDUCTION:\n")
cat("     Compared with treat-all (scan everyone), using the nomogram\n")
cat("     at a threshold where net benefit is superior implies fewer\n")
cat("     false-positive referrals per true metastasis detected.\n\n")

cat("  4. LIMITATIONS:\n")
cat("     - DCA is based on apparent performance (same dataset used\n")
cat("       for model development); external validation is required\n")
cat("       before clinical deployment.\n")
cat("     - Net benefit estimates at extreme thresholds are unstable\n")
cat("       with n = 114.\n")
cat("     - Clinical acceptability of threshold range requires\n")
cat("       clinician input on acceptable miss rate vs. scan burden.\n")


# =============================================================================
# SECTION 6 — PART 7: Publication-Ready Summary Paragraph
# =============================================================================

cat("\n=== PART 7: PUBLICATION-READY PARAGRAPH ===\n\n")

if (!is.na(threshold_lower)) {
  cat(sprintf(paste0(
    "Decision curve analysis was performed to evaluate the clinical utility ",
    "of the prediction nomogram for occult distant metastasis in gastric cancer ",
    "across threshold probabilities ranging from 5%% to 50%%. ",
    "The nomogram demonstrated greater net benefit than both the treat-all and ",
    "treat-none strategies across threshold probabilities of approximately ",
    "%.0f%% to %.0f%%, indicating potential clinical utility for prioritizing ",
    "advanced imaging (CT/MRI) in this probability range. ",
    "Within this threshold interval, use of the nomogram may reduce unnecessary ",
    "CT/MRI referrals compared with a strategy of universal imaging, while ",
    "maintaining superior identification of patients at elevated risk of occult ",
    "distant metastasis compared with a no-imaging strategy. ",
    "These findings support the nomogram as a decision support tool for ",
    "clinicians determining which patients evaluated by endoscopy, ultrasound, ",
    "and complete blood count may benefit from prioritized advanced imaging. ",
    "External validation is required before clinical deployment."
  ),
  threshold_lower * 100,
  threshold_upper * 100
  ))
} else {
  cat(paste0(
    "Decision curve analysis was performed across threshold probabilities ",
    "from 5% to 50%. The nomogram did not demonstrate consistent net benefit ",
    "superior to default strategies across this range. These findings suggest ",
    "limited clinical utility of the current model for CT/MRI prioritization ",
    "and indicate that further model refinement or external validation is ",
    "required before clinical application."
  ))
}

cat("\n")


# =============================================================================
# SECTION 7: Step 4 Summary
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("STEP 4 COMPLETE — DECISION CURVE ANALYSIS SUMMARY\n")
cat("================================================================\n")
cat(sprintf("Model:       Linear logistic regression (approved Steps 1-3)\n"))
cat(sprintf("n = %d | Events = %d | Prevalence = %.1f%%\n",
            nrow(df), sum(df$Metastasis_TNM),
            mean(df$Metastasis_TNM) * 100))
cat(sprintf("Threshold range evaluated: 5%% to 50%%\n"))
if (!is.na(threshold_lower)) {
  cat(sprintf("Net benefit superior to defaults: ~%.0f%% to %.0f%%\n",
              threshold_lower * 100, threshold_upper * 100))
} else {
  cat("Net benefit: nomogram did not consistently outperform defaults\n")
}
cat("---\n")
cat("Figures:  figure6_predicted_prob_distribution.png\n")
cat("          figure7_DCA.png\n")
cat("Tables:   Table 5 (net benefit by threshold)\n")
cat("================================================================\n")
cat("STOPPING — Awaiting review before Step 5 (Elastic Net Benchmark).\n")