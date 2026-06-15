#!/usr/bin/env Rscript

############################################################
## Required packages
############################################################
# tidymodels:
#   modeling workflow, recipes, resampling, and metrics
#
# ranger:
#   random forest engine used through tidymodels
#
# dplyr:
#   data manipulation, filtering, mutation, and column selection
library(tidymodels)
library(ranger)
library(dplyr)


############################################################
## Command-line arguments
############################################################
# Usage:
#   Rscript rf_cv_allstats.r \
#       <input_file> \
#       <output_dir> \
#       <dataset_name> \
#       <evol_mode_comparison> \
#       <positive_class> \
#       <negative_class> \
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
#   Rscript rf_cv_allstats.r \
#       moments_over_vs_neutral.txt \
#       results_rf \
#       moments \
#       over_vs_neutral \
#       Overdominance \
#       Neutral \
#       2Nes=500 \
#       results_all_stats_rf_cv10_over_vs_neutral_moments.txt

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 8) {
  stop(
    "Usage: Rscript rf_cv_allstats.r ",
    "<input_file> <output_dir> <dataset_name> <evol_mode_comparison> ",
    "<positive_class> <negative_class> <selection_strength> ",
    "<metrics_output_file>"
  )
}

input_file <- args[1]
out_dir <- args[2]
dataset_name <- args[3]
evol_mode_comparison <- args[4]
positive_class <- args[5]
negative_class <- args[6]
selection_strength <- args[7]
metrics_output_file <- args[8]


############################################################
## Create output directories
############################################################
# Main output directory.
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Directory for fold-level RF predictions.
prediction_dir <- file.path(out_dir, "rf_predictions_all_stats")
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
## Function: RF classification using all statistics
############################################################
# This function runs one RF analysis for:
#   - one input dataset
#   - one selection-strength value
#   - one evolutionary comparison
#   - all available summary statistics together
#
# Required metadata columns:
#   Selection_mode
#   Selection_strength
#   Run
#
# Predictor columns:
#   all remaining columns after removing metadata columns
#
# Cross-validation:
#   - grouped 10-fold CV
#   - grouping variable = Run
#   - prevents samples from the same simulation replicate
#     from being split across training/test folds
#
# Positive/first class:
#   controlled by the positive_class and negative_class arguments
#
# Outputs:
#   1. fold-level prediction file
#   2. cross-validation metrics returned as a tibble

run_all_stats_cv <- function(df,
                             selection_strength,
                             dataset_name,
                             evol_mode_comparison,
                             positive_class,
                             negative_class,
                             out_dir) {

  message(
    "Running RF all-statistics analysis: ",
    "Selection_strength = ", selection_strength,
    " | dataset = ", dataset_name,
    " | comparison = ", evol_mode_comparison
  )


  ############################################################
  ## Select predictor columns
  ############################################################
  # All columns except metadata columns are used as predictors.

  metadata_cols <- c(
    "Selection_mode",
    "Selection_strength",
    "Run"
  )

  missing_metadata <- setdiff(metadata_cols, colnames(df))

  if (length(missing_metadata) > 0) {
    stop(
      "Input file is missing required metadata column(s): ",
      paste(missing_metadata, collapse = ", ")
    )
  }

  stat_l <- setdiff(colnames(df), metadata_cols)

  if (length(stat_l) == 0) {
    stop("No predictor columns found after removing metadata columns.")
  }


  ############################################################
  ## Subset data
  ############################################################
  # Keep only one selection-strength category.
  # Run is retained temporarily for grouped cross-validation.

  df_sub <- df %>%
    filter(Selection_strength == selection_strength) %>%
    dplyr::select(
      Selection_mode,
      Selection_strength,
      Run,
      all_of(stat_l)
    )

  if (nrow(df_sub) == 0) {
    stop(
      "No rows found for Selection_strength = ",
      selection_strength
    )
  }

  # Define class order explicitly.
  # The first factor level is treated as the event level by yardstick
  # unless event_level is otherwise specified.
  df_sub$Selection_mode <- factor(
    df_sub$Selection_mode,
    levels = c(positive_class, negative_class)
  )


  ############################################################
  ## Grouped 10-fold cross-validation
  ############################################################
  # Grouping by Run prevents data from the same simulation replicate
  # from being split across training and test folds.

  set.seed(123)

  folds <- group_vfold_cv(
    df_sub,
    v = 10,
    group = Run
  )

  # Save predictions and workflow objects from each fold.
  ctrl <- control_resamples(
    save_pred = TRUE,
    save_workflow = TRUE
  )


  ############################################################
  ## Remove metadata columns before modeling
  ############################################################
  # Selection_strength and Run should not be used as predictors.

  df_model <- df_sub %>%
    dplyr::select(-Selection_strength, -Run)


  ############################################################
  ## Preprocessing recipe
  ############################################################
  # step_zv():
  #   remove predictors with zero variance
  #
  # step_normalize():
  #   center and scale predictors

  rec <- recipe(Selection_mode ~ ., data = df_model) %>%
    step_zv(all_predictors()) %>%
    step_normalize(all_predictors())


  ############################################################
  ## Random forest model
  ############################################################
  # trees:
  #   number of trees
  #
  # mtry:
  #   number of predictors randomly sampled at each split
  #   here set to floor(sqrt(number of predictors))
  #
  # min_n:
  #   minimum node size

  rf_mod <- rand_forest(
    mode = "classification",
    trees = 500,
    mtry = floor(sqrt(length(stat_l))),
    min_n = 5
  ) %>%
    set_engine("ranger")


  ############################################################
  ## Workflow
  ############################################################
  # Combines preprocessing recipe and RF model.

  wf <- workflow() %>%
    add_model(rf_mod) %>%
    add_recipe(rec)


  ############################################################
  ## Fit model using grouped CV
  ############################################################

  cv_results <- wf %>%
    fit_resamples(
      resamples = folds,
      metrics = metric_set(
        accuracy,
        pr_auc,
        precision,
        recall,
        roc_auc
      ),
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
      statistic = "all_stats"
    )

  prediction_output_file <- file.path(
    out_dir,
    "rf_predictions_all_stats",
    paste0(
      "rf_predictions_all_stats_",
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
  ## Collect and return CV summary metrics
  ############################################################

  summarised <- cv_results %>%
    collect_metrics() %>%
    mutate(
      statistic = "all_stats",
      Selection_strength = selection_strength,
      dataset = dataset_name,
      evol_mode_comparison = evol_mode_comparison
    )

  message("Finished RF all-statistics analysis.")
  message("Predictions written to: ", prediction_output_file)

  return(summarised)
}


############################################################
## Run analysis
############################################################

results <- run_all_stats_cv(
  df = df,
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