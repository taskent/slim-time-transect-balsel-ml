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
## Output directory
############################################################
out_dir <- "output_directory"  # update with actual path


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
#   - factor level order is:
#       Overdominance, Neutral (Overdominance, NFDS for the other comparison)
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

  # Define class order explicitly.
  df_sub$Selection_mode <- factor(
    df_sub$Selection_mode,
    levels = c("Overdominance", "Neutral")
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
    roc_auc
    pr_auc,
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

  write.table(
    preds,
    file = file.path(
      out_dir,
      paste0(
        "svm_predictions_ind_stats/svm_predictions_",
        summ_stat_name, "_",
        dataset_name, "_",
        selection_strength, "_",
        evol_mode_comparison,
        ".txt"
      )
    ),
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

  return(summarised)
}