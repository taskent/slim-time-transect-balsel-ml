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
############################################################
## Command-line arguments
############################################################
# Usage:
#   Rscript rf_cv_singlestat.r \
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
#   Rscript rf_cv_singlestat.r \
#       moments_over_vs_neutral.txt \
#       results_rf \
#       moments \
#       over_vs_neutral \
#       Overdominance \
#       Neutral \
#       TjD \
#       2Nes=500 \
#       results_single_stats_rf_cv10_over_vs_neutral_moments.txt

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 9) {
  stop(
    "Usage: Rscript rf_cv_singlestat.r ",
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
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

prediction_dir <- file.path(out_dir, "rf_predictions_ind_stats")
dir.create(prediction_dir, showWarnings = FALSE, recursive = TRUE)


############################################################
## Read input dataset
############################################################
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
#   base summary-statistic name, e.g. TjD, Median_R2, Max_B0maf
#
# dataset_name:
#   last_sampling:
#       uses one predictor: summ_stat_name
#
#   moments:
#       uses four predictors:
#       mean, variance, skewness, kurtosis
#
#   autocorrelation:
#       uses five predictors:
#       autocorrelation lags 1 to 5
#
#   moments_autocorrelation:
#       uses nine predictors:
#       four moment predictors + five autocorrelation predictors

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
## Function: RF classification for one statistic
############################################################
run_single_stat_cv <- function(df,
                               summ_stat_name,
                               selection_strength,
                               dataset_name,
                               evol_mode_comparison,
                               positive_class,
                               negative_class,
                               out_dir) {

  message(
    "Running RF single-stat: ", summ_stat_name,
    " | Selection_strength = ", selection_strength,
    " | dataset = ", dataset_name,
    " | evol_mode_comparison = ", evol_mode_comparison
  )

  ############################################################
  ## Determine predictor columns
  ############################################################
  
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
  # Keep only one selection-strength class.
  # Run is retained temporarily for grouped cross-validation.

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

  ctrl <- control_resamples(
    save_pred = TRUE,
    save_workflow = TRUE
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
  # For one predictor, mtry = 1.
  # For multiple derived predictors, mtry = floor(sqrt(p)).

  rf_mod <- rand_forest(
    mode = "classification",
    trees = 500,
    mtry = floor(sqrt(length(new_stat_name_l))),
    min_n = 5
  ) %>%
    set_engine("ranger")

  ############################################################
  ## Workflow
  ############################################################

  wf <- workflow() %>%
    add_model(rf_mod) %>%
    add_recipe(rec)

  ############################################################
  ## Fit model with grouped CV
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
      statistic = summ_stat_name,
      Selection_strength = selection_strength,
      dataset = dataset_name,
      evol_mode_comparison = evol_mode_comparison
    )

  prediction_output_file <- file.path(
    out_dir,
    "rf_predictions_ind_stats",
    paste0(
      "rf_predictions_",
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
  ## Collect CV summary metrics
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
    "Script finished: ", summ_stat_name,
    " | Selection_strength = ", selection_strength,
    " | dataset = ", dataset_name,
    " | evol_mode_comparison = ", evol_mode_comparison
  )
  message("Predictions written to: ", prediction_output_file)

  return(summarised)
}


############################################################
## Run analysis
############################################################
results <- run_single_stat_cv(
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


