############################################################
## Required packages
############################################################
# tidymodels: modeling workflow, resampling, metrics
# ranger: random forest engine
# dplyr: data manipulation
library(tidymodels)
library(ranger)
library(dplyr)


############################################################
## Output directory
############################################################
out_dir <- "output_directory" # Update with actual path


############################################################
## Function: RF classification using all statistics
############################################################
# This function runs a random forest classifier for:
#   - one input dataset (i.e., either Overdominance vs. Neutral or Overdominance vs. NFDS; and for each of these two comparisons, either (1) last sampling, (2) moments, (3) autocorrelation or (4) moments + autocorrelation data sets)
#   - one selection-strength value (e.g., 2Nes=100)
#   - all selected summary statistics together
#
# Cross-validation:
#   - 10-fold grouped CV
#   - grouping variable = Run
#
# Positive/first class:
#   - factor level order is:
#       Overdominance, Neutral (Overdominance, NFDS for the other comparison)
#
# Outputs:
#   1. Prediction file for each Selection_strength × dataset
#   2. Summary metrics returned as a tibble

run_all_stats_cv <- function(df, stat_l, selection_strength, dataset_name,
                             evol_mode_comparison, ind_or_all_stats) {
  
  message("Running: Selection_strength = ", selection_strength,
          " | dataset = ", dataset_name,
          " | evol_mode_comparison = ", evol_mode_comparison,
          " | ind_or_all_stats = all statistics together")
  
  ############################################################
  ## Subset data
  ############################################################
  # Keep only one Selection_strength value and selected predictor columns.
  # Run is retained temporarily for grouped cross-validation.
  df_sub <- df %>% 
    filter(Selection_strength == selection_strength) %>% 
    dplyr::select(Selection_mode, Selection_strength, Run, all_of(stat_l))

  # Define class order explicitly.
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
  
  # Save predictions and workflow objects from each fold.
  ctrl <- control_resamples(
    save_pred = TRUE,
    save_workflow = TRUE
  )

  ############################################################
  ## Remove metadata columns before modeling
  ############################################################
  # Selection_strength and Run should not be used as predictors.
  drops <- c("Selection_strength", "Run")
  df_sub <- df_sub[, !(names(df_sub) %in% drops)]

  print(colnames(df_sub))
  
  ############################################################
  ## Recipe
  ############################################################
  # step_zv(): remove predictors with zero variance.
  # step_normalize(): standardize predictors.
  rec <- recipe(Selection_mode ~ ., data = df_sub) %>%
    step_zv(all_predictors()) %>%
    step_normalize(all_predictors())
  
  ############################################################
  ## Random forest model
  ############################################################
  # trees = 500
  # mtry = sqrt(number of predictors), rounded down
  # min_n = minimum node size
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
  wf <- workflow() %>%
    add_model(rf_mod) %>%
    add_recipe(rec)
  
  ############################################################
  ## Fit model using grouped CV
  ############################################################
  cv_results <- wf %>% fit_resamples(
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
      ind_or_all_stats = ind_or_all_stats
    )
    
  write.table(
    preds,
    file = file.path(
      out_dir,
      paste0(
        "rf_predictions_allstats_",
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
  ## Collect and return CV metrics
  ############################################################
  summarised <- cv_results %>% 
    collect_metrics() %>%
    mutate(
      statistic = "all_stats",
      Selection_strength = selection_strength,
      dataset = dataset_name,
      evol_mode_comparison = evol_mode_comparison
    )

  message("Script finished: Selection_strength = ", selection_strength,
          " | dataset = ", dataset_name,
          " | evol_mode_comparison = ", evol_mode_comparison,
          " | ind_or_all_stats = all statistics together")
  
  return(summarised)
}