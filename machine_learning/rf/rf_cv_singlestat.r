############################################################
## Required packages
############################################################
library(tidymodels)  # workflows, recipes, resampling, metrics
library(ranger)      # random forest engine
library(dplyr)       # data manipulation


############################################################
## Output directory
############################################################
out_dir <- "output_directory"  # update with actual path


############################################################
## Helper function: choose predictor columns by dataset type
############################################################
# dataset_name should be one of:
#   "last_sampling"
#   "moments"
#   "autocorrelation"
#   "moments_autocorrelation"
#
# summ_stat_name is the base name of the statistic, e.g. "TjD".
#
# This returns the actual column name(s) used as RF predictors.

get_predictor_columns <- function(summ_stat_name, dataset_name) {

  if (dataset_name == "last_sampling") {

    new_stat_name_l <- summ_stat_name

  } else if (dataset_name == "moments") {

    new_stat_name_l <- c(
      paste0(summ_stat_name, "_mean"),
      paste0(summ_stat_name, "_sd"),
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
      paste0(summ_stat_name, "_sd"),
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
# Runs RF classification with grouped 10-fold CV for:
#   - one summary statistic
#   - one selection-strength value
#   - one dataset type
#
# The actual predictor columns are determined by dataset_name:
#   last_sampling             -> 1 predictor
#   moments                   -> 4 predictors
#   autocorrelation           -> 5 predictors
#   moments_autocorrelation   -> 9 predictors
#
# Positive/first class:
#   - factor level order is:
#       Overdominance, Neutral (Overdominance, NFDS for the other comparison)
#
# Grouped CV:
#   group = Run
#   This prevents observations from the same simulation replicate
#   from being split across training and test folds.

run_single_stat_cv <- function(df,
                               summ_stat_name,
                               selection_strength,
                               dataset_name,
                               evol_mode_comparison) {

  message(
    "Running: ", summ_stat_name,
    " | Selection_strength = ", selection_strength,
    " | dataset = ", dataset_name,
    " | evol_mode_comparison = ", evol_mode_comparison
  )

  ############################################################
  ## Determine predictor columns for this dataset type
  ############################################################
  new_stat_name_l <- get_predictor_columns(
    summ_stat_name = summ_stat_name,
    dataset_name = dataset_name
  )

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

  # Define class order explicitly.
  # The first level is treated as the event level by yardstick
  # unless event_level is otherwise specified.
  df_sub$Selection_mode <- factor(
    df_sub$Selection_mode,
    levels = c("Overdominance", "Neutral")
  )

  ############################################################
  ## Grouped 10-fold cross-validation
  ############################################################
  set.seed(123)

  folds <- group_vfold_cv(
    df_sub,
    v = 10,
    group = Run
  )

  # Save fold-level predictions and workflows.
  ctrl <- control_resamples(
    save_pred = TRUE,
    save_workflow = TRUE
  )

  ############################################################
  ## Remove metadata columns before model fitting
  ############################################################
  # Selection_strength and Run should not be used as predictors.
  df_sub <- df_sub %>%
    dplyr::select(-Selection_strength, -Run)

  ############################################################
  ## Preprocessing recipe
  ############################################################
  # step_zv(): remove zero-variance predictors
  # step_normalize(): center and scale predictors
  rec <- recipe(Selection_mode ~ ., data = df_sub) %>%
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

  write.table(
    preds,
    file = file.path(
      out_dir,
      paste0(
        "rf_predictions_ind_stats/rf_predictions_",
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
  ## Return CV summary metrics
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

  return(summarised)
}
