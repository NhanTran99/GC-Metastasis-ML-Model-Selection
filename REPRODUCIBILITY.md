# Reproducibility Guide

## Software

- R version 4.4.1 
- rms
- pROC
- rmda
- glmnet
- xgboost
- tidyverse

## Workflow

Run scripts sequentially:

1. `01_logistic_regression.R`
2. `02_bootstrap_validation.R`
3. `03_rcs_analysis.R`
4. `04_decision_curve_analysis.R`
5. `05_elastic_net_benchmark.R`
6. `06_xgboost_benchmark.R`
7. `07_final_model_selection.R`

All analyses were conducted in R.

## Data

Original dataset is not publicly available due to confidentiality restrictions.

Users may adapt the scripts to compatible datasets.