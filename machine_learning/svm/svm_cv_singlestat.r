#!/usr/bin/env Rscript

############################################################
## Required packages
############################################################
# tidymodels: modeling framework, resampling, metrics
# kernlab: SVM backend
# dplyr: data manipulation
# dials: parameter-related infrastructure used by tidymodels
# recipes: preprocessing
# workflows: combine recipe + model
# yardstick: classification metrics
library(tidymodels)
library(kernlab)
library(dplyr)
library(dials)

library(recipes)
library(workflows)
library(yardstick)

options(tidymodels.dark = TRUE)


############################################################
## Command-line arguments
############################################################
# Usage:
#   Rscript svm_cv_singlestat.r \
#       <input_file> \
#       <output_dir> \
#       <dataset_name> \
#       <evol_mode_comparison> \
#       <positive_class> \
#       <negative_class> \
#       <summ_stat_name> \
#       <selection_strength> \
#       <metrics_output_file>
#
# dataset_name should be one of:
#   last_sampling
#   moments
#   autocorrelation
#   moments_autocorrelation
#
# Example:
#   Rscript svm_cv_singlestat.r \
#       moments_over_vs_neutral.txt \
#       results_svm \
#       moments \
#       over_vs_neutral \
#       Overdominance \
#       Neutral \
#       TjD \
#       2Nes=500 \
#       results_single_stats_svm_cv10_over_vs_neutral_moments.txt

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 9) {
  stop(
    "Usage: Rscript svm_cv_singlestat.r ",
    "<input_file> <output_dir> <dataset_name> <evol_mode_comparison> ",
    "<positive_class> <negative_class> <summ_stat_name> ",
    "<selection_strength> <metrics_output_file>"
  )
}

input_file <- args[1]
out_dir <- args[2]
dataset_name <- args[3]
evol_mode_comparison <- args[4]
positive_class <- args[5]
negative_class <- args[6]
summ_stat_name <- args[7]
selection_strength <- args[8]
metrics_output_file <- args[9]


############################################################
## Create output directories
############################################################
# Main output directory.
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Directory for fold-level SVM predictions.
prediction_dir <- file.path(out_dir, "svm_predictions_ind_stats")
dir.create(prediction_dir, showWarnings = FALSE, recursive = TRUE)


############################################################
## Read input dataset
############################################################
# Input file should contain:
#   Selection_mode
#   Selection_strength
#   Run
#   summary-statistic predictor columns

df <- read.table(
  input_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)


############################################################
## Helper function: choose predictor columns by dataset type
############################################################
# summ_stat_name:
#   base summary-statistic name, e.g. "TjD", "Median_R2", "Max_B0maf"
#
# dataset_name:
#   "last_sampling"
#   "moments"
#   "autocorrelation"
#   "moments_autocorrelation"
#
# Output:
#   new_stat_name_l: vector of actual predictor-column names
#
# Examples:
#   last_sampling:
#       "TjD"
#
#   moments:
#       "TjD_mean", "TjD_var", "TjD_skewness", "TjD_kurtosis"
#
#   autocorrelation:
#       "autocorrLag1_TjD", ..., "autocorrLag5_TjD"
#
#   moments_autocorrelation:
#       all four moments + all five autocorrelation features

get_predictor_columns <- function(summ_stat_name, dataset_name) {

  if (dataset_name == "last_sampling") {

    new_stat_name_l <- summ_stat_name

  } else if (dataset_name == "moments") {

    new_stat_name_l <- c(
      paste0(summ_stat_name, "_mean"),
      paste0(summ_stat_name, "_var"),
      paste0(summ_stat_name, "_skewness"),
      paste0(summ_stat_name, "_kurtosis")
    )

  } else if (dataset_name == "autocorrelation") {

    new_stat_name_l <- c(
      paste0("autocorrLag1_", summ_stat_name),
      paste0("autocorrLag2_", summ_stat_name),
      paste0("autocorrLag3_", summ_stat_name),
      paste0("autocorrLag4_", summ_stat_name),
      paste0("autocorrLag5_", summ_stat_name)
    )

  } else if (dataset_name == "moments_autocorrelation") {

    new_stat_name_l <- c(
      paste0(summ_stat_name, "_mean"),
      paste0(summ_stat_name, "_var"),
      paste0(summ_stat_name, "_skewness"),
      paste0(summ_stat_name, "_kurtosis"),
      paste0("autocorrLag1_", summ_stat_name),
      paste0("autocorrLag2_", summ_stat_name),
      paste0("autocorrLag3_", summ_stat_name),
      paste0("autocorrLag4_", summ_stat_name),
      paste0("autocorrLag5_", summ_stat_name)
    )

  } else {
    stop(
      "dataset_name must be one of: last_sampling, moments, ",
      "autocorrelation, moments_autocorrelation"
    )
  }

  return(new_stat_name_l)
}


############################################################
## Function: linear SVM classification for one statistic
############################################################
# This function runs a linear SVM classifier for:
#   - one input data frame
#   - one base summary statistic
#   - one selection-strength value
#   - one dataset type
#
# The dataset type determines which predictor columns are used:
#   last_sampling             -> 1 predictor
#   moments                   -> 4 predictors
#   autocorrelation           -> 5 predictors
#   moments_autocorrelation   -> 9 predictors
#
# Cross-validation:
#   - grouped 10-fold CV
#   - grouping variable = Run
#   - this prevents samples from the same simulation replicate
#     being split across train/test folds.
#
# Positive/first class:
#   - controlled by positive_class and negative_class arguments
#
# Class weights:
#   - computed from class frequencies within the current subset
#   - passed to kernlab to reduce class-imbalance effects.
#
# Outputs:
#   1. fold-level predictions written to file
#   2. cross-validation metrics returned as a tibble

run_svm_single_stat_cv <- function(df,
                                   summ_stat_name,
                                   selection_strength,
                                   dataset_name,
                                   evol_mode_comparison,
                                   positive_class,
                                   negative_class,
                                   out_dir) {

  message("--------------------------------------------------")
  message(
    "Running SVM: dataset = ", dataset_name,
    " | statistic = ", summ_stat_name,
    " | Selection_strength = ", selection_strength
  )
  message("--------------------------------------------------")


  ############################################################
  ## Determine predictor columns
  ############################################################
  # This replaces the separate last_sampling / moments /
  # autocorrelation / moments_autocorrelation functions.
  new_stat_name_l <- get_predictor_columns(
    summ_stat_name = summ_stat_name,
    dataset_name = dataset_name
  )

  missing_cols <- setdiff(new_stat_name_l, colnames(df))

  if (length(missing_cols) > 0) {
    stop(
      "Missing predictor column(s): ",
      paste(missing_cols, collapse = ", ")
    )
  }


  ############################################################
  ## Subset data
  ############################################################
  # Keep only one selection-strength value.
  # Run is kept temporarily for grouped CV.
  df_sub <- df %>%
    filter(Selection_strength == selection_strength) %>%
    dplyr::select(
      Selection_mode,
      Selection_strength,
      Run,
      all_of(new_stat_name_l)
    )

  if (nrow(df_sub) == 0) {
    stop("No rows found for Selection_strength = ", selection_strength)
  }

  ############################################################
  ## Remove rows with non-finite predictor values
  ############################################################
  # Rows containing NA, NaN, Inf or -Inf in any predictor are removed.
  # This avoids prediction-row mismatches during resampling.
  n_before_filter <- nrow(df_sub)

  df_sub <- df_sub %>%
    dplyr::filter(
      if_all(
        -c(Selection_mode, Selection_strength, Run),
        ~ is.finite(.x)
      )
    )

  n_after_filter <- nrow(df_sub)

  message(
    "Removed ",
    n_before_filter - n_after_filter,
    " rows with non-finite predictor values."
  )

  if (n_after_filter == 0) {
    stop(
      "No rows remain after removing non-finite predictor values for Selection_strength = ",
      selection_strength
    )
  }

  # Define class order explicitly.
  df_sub$Selection_mode <- factor(
    df_sub$Selection_mode,
    levels = c(positive_class, negative_class)
  )


  ############################################################
  ## Class weights
  ############################################################
  # Weight = total sample size / (number of classes × class count)
  # This upweights the minority class.
  class_counts <- table(df_sub$Selection_mode)
  weights <- sum(class_counts) / (2 * class_counts)


  ############################################################
  ## Grouped 10-fold cross-validation
  ############################################################
  set.seed(123)

  folds <- group_vfold_cv(
    df_sub,
    v = 10,
    group = Run
  )


  ############################################################
  ## Remove metadata columns before modeling
  ############################################################
  # Selection_strength and Run are not predictors.
  df_model <- df_sub %>%
    dplyr::select(-Selection_strength, -Run)


  ############################################################
  ## Preprocessing recipe
  ############################################################
  # step_zv(): removes zero-variance predictors
  # step_normalize(): centers and scales predictors

  rec <- recipe(Selection_mode ~ ., data = df_model) %>%
    step_zv(all_predictors()) %>%
    step_normalize(all_predictors())


  ############################################################
  ## Linear SVM model
  ############################################################
  # This script uses a fixed linear SVM with no hyperparameter tuning.
  # Class weights are passed to kernlab.
  svm_mod <- svm_linear() %>%
    set_engine("kernlab", class.weights = weights) %>%
    set_mode("classification")


  ############################################################
  ## Workflow
  ############################################################
  # Combines preprocessing recipe and SVM model.
  wf <- workflow() %>%
    add_recipe(rec) %>%
    add_model(svm_mod)


  ############################################################
  ## Metrics
  ############################################################
  # Metrics are calculated for each cross-validation fold.
  metrics_used <- metric_set(
    accuracy,
    precision,
    recall,
    roc_auc,
    pr_auc
  )


  ############################################################
  ## Resampling control
  ############################################################
  # save_pred = TRUE:
  #   keep fold-level predictions for later inspection.
  #
  # save_workflow = TRUE:
  #   store workflow information.
  ctrl <- control_resamples(
    save_pred = TRUE,
    save_workflow = TRUE
  )


  ############################################################
  ## Fit model with grouped CV
  ############################################################
  cv_results <- wf %>%
    fit_resamples(
      resamples = folds,
      metrics = metrics_used,
      control = ctrl
    )


  ############################################################
  ## Save fold-level predictions
  ############################################################
  preds <- cv_results %>%
    collect_predictions() %>%
    mutate(
      Selection_strength = selection_strength,
      dataset = dataset_name,
      evol_mode_comparison = evol_mode_comparison,
      statistic = summ_stat_name
    )

  prediction_output_file <- file.path(
    out_dir,
    "svm_predictions_ind_stats",
    paste0(
      "svm_predictions_",
      summ_stat_name, "_",
      dataset_name, "_",
      selection_strength, "_",
      evol_mode_comparison,
      ".txt"
    )
  )

  write.table(
    preds,
    file = prediction_output_file,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )


  ############################################################
  ## Collect and return summary metrics
  ############################################################
  summarised <- cv_results %>%
    collect_metrics() %>%
    mutate(
      statistic = summ_stat_name,
      Selection_strength = selection_strength,
      dataset = dataset_name,
      evol_mode_comparison = evol_mode_comparison
    )

  message(
    "Finished SVM run: evol_mode_comparison = ", evol_mode_comparison,
    " | dataset = ", dataset_name,
    " | statistic = ", summ_stat_name,
    " | Selection_strength = ", selection_strength
  )

  message("Predictions written to: ", prediction_output_file)

  return(summarised)
}


############################################################
## Run analysis
############################################################

results <- run_svm_single_stat_cv(
  df = df,
  summ_stat_name = summ_stat_name,
  selection_strength = selection_strength,
  dataset_name = dataset_name,
  evol_mode_comparison = evol_mode_comparison,
  positive_class = positive_class,
  negative_class = negative_class,
  out_dir = out_dir
)


############################################################
## Append metrics to output file
############################################################
# If the output file does not exist yet, write column names.
# If it already exists, append rows without column names.

write_colnames <- !file.exists(metrics_output_file)

write.table(
  results,
  file = metrics_output_file,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE,
  col.names = write_colnames,
  append = !write_colnames
)


############################################################
## Finished
############################################################

message("Metrics appended to: ", metrics_output_file)