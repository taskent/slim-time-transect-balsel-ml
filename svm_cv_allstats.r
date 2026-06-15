############################################################
## Required packages
############################################################
# tidymodels: model specifications, workflows, resampling, metrics
# kernlab: SVM backend used by tidymodels
# dplyr: data manipulation
# dials: tuning parameter definitions and parameter ranges
# recipes: preprocessing recipe steps
# workflows: model + recipe workflows
# tune: grid tuning
# yardstick: performance metrics
library(tidymodels)
library(kernlab)
library(dplyr)
library(dials)

library(recipes)
library(workflows)
library(tune)
library(yardstick)

options(tidymodels.dark = TRUE)


############################################################
## Output directory
############################################################
out_dir <- "output_directory"  # update with actual path


############################################################
## Function: SVM grid tuning with grouped 10-fold CV
############################################################
# This function runs SVM classification for:
#   - one dataset type
#   - one selection-strength value
#   - one kernel type
#   - a chosen set of summary statistics
#
# Supported kernels:
#   - linear
#   - polynomial
#   - radial
#
# Cross-validation:
#   - 10-fold grouped CV
#   - grouping variable = Run
#
# Positive/first class:
#   - factor level order is:
#       Overdominance, Neutral (Overdominance, NFDS for the other comparison)
#
# Class weights:
#   - calculated from class frequencies within the current
#     selection-strength subset
#   - passed to kernlab to reduce effects of class imbalance
#
# Tuning:
#   - the full grid is evaluated using multiple metrics
#   - the best hyperparameter combination is selected using
#     ROC AUC
#
# Outputs:
#   1. Full tuning object saved as RDS
#   2. Tuning notes saved as TXT
#   3. Best hyperparameters saved as TXT
#   4. Metrics table returned as a tibble

run_svm_grid_cv <- function(df,
                            evol_mode_comparision,
                            dataset_name,
                            stat_l,
                            selection_strength,
                            kernel_mode,
                            out_dir) {

  message("--------------------------------------------------")
  message(Sys.time(), " Starting SVM tuning")
  message(
    "Running SVM grid tuning: Selection_strength = ", selection_strength,
    " | kernel = ", kernel_mode,
    " | dataset = ", dataset_name
  )
  message("--------------------------------------------------")


  ############################################################
  ## Subset data
  ############################################################
  # Keep only the requested selection-strength value.
  # Run is retained temporarily for grouped cross-validation.
  df_sub <- df %>%
    filter(Selection_strength == selection_strength) %>%
    dplyr::select(
      Selection_mode,
      Selection_strength,
      Run,
      all_of(stat_l)
    )

  # Define class order explicitly.
  df_sub$Selection_mode <- factor(
    df_sub$Selection_mode,
    levels = c("Overdominance", "Neutral")
  )


  ############################################################
  ## Class weights
  ############################################################
  # Class weights are inversely proportional to class counts.
  # This gives larger weight to the minority class.
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
  # Selection_strength and Run should not be used as predictors.
  df_model <- df_sub %>%
    dplyr::select(-Selection_strength, -Run)


  ############################################################
  ## Preprocessing recipe
  ############################################################
  # step_zv(): removes zero-variance predictors
  # step_normalize(): centers and scales predictors
  #
  # Normalization is important for SVMs because SVM kernels are
  # distance-based / inner-product-based.
  rec <- recipe(Selection_mode ~ ., data = df_model) %>%
    step_zv(all_predictors()) %>%
    step_normalize(all_predictors())


  ############################################################
  ## Model specification and tuning parameter space
  ############################################################

  if (kernel_mode == "linear") {

    # Linear SVM:
    # tuning parameter = cost
    svm_mod <- svm_linear(
      cost = tune()
    ) %>%
      set_engine("kernlab", class.weights = weights) %>%
      set_mode("classification")

    params <- parameters(
      cost(range = c(-5, 5))
    )

  } else if (kernel_mode == "polynomial") {

    # Polynomial SVM:
    # tuning parameters:
    #   cost
    #   polynomial degree
    #   scale factor
    svm_mod <- svm_poly(
      cost = tune(),
      degree = tune("prod_degree"),
      scale_factor = tune()
    ) %>%
      set_engine("kernlab", class.weights = weights) %>%
      set_mode("classification")

    params <- parameters(
      cost(range = c(-5, 5)),
      scale_factor(range = c(-3, -1)),
      prod_degree(range = c(1L, 3L))
    )

  } else if (kernel_mode == "radial") {

    # Radial basis function SVM:
    # tuning parameters:
    #   cost
    #   rbf_sigma
    svm_mod <- svm_rbf(
      cost = tune(),
      rbf_sigma = tune()
    ) %>%
      set_engine("kernlab", class.weights = weights) %>%
      set_mode("classification")

    params <- parameters(
      cost(range = c(-5, 5)),
      rbf_sigma(range = c(-3, -1))
    )

  } else {
    stop("kernel_mode must be one of: linear, polynomial, radial")
  }


  ############################################################
  ## Workflow
  ############################################################
  # Combines preprocessing recipe and model specification.
  wf <- workflow() %>%
    add_recipe(rec) %>%
    add_model(svm_mod)


  ############################################################
  ## Metrics
  ############################################################
  # Metrics saved for each parameter combination and CV fold.
  #
  # These metrics are all calculated, but the final "best"
  # hyperparameter combination is selected later using ROC AUC.
  metrics_used <- metric_set(
    accuracy, 
    precision,
    recall,
    pr_auc,
    roc_auc
  )


  ############################################################
  ## Tuning control
  ############################################################
  # save_pred = TRUE: save fold-level predictions
  # save_workflow = TRUE: save workflow information
  #
  # Parallel settings are intentionally removed.
  ctrl <- control_grid(
    verbose = TRUE,
    save_pred = TRUE,
    save_workflow = TRUE
  )


  ############################################################
  ## Build tuning grid
  ############################################################
  # Grid size differs by kernel complexity.
  if (kernel_mode == "linear") {
    grid_vals <- grid_space_filling(params, size = 20)
  } else if (kernel_mode == "polynomial") {
    grid_vals <- grid_space_filling(params, size = 60)
  } else if (kernel_mode == "radial") {
    grid_vals <- grid_space_filling(params, size = 40)
  }


  ############################################################
  ## Run grid tuning
  ############################################################
  message(Sys.time(), " Starting tune_grid()")

  svm_tuned <- tryCatch(
    {
      tune_grid(
        wf,
        resamples = folds,
        grid = grid_vals,
        metrics = metrics_used,
        control = ctrl
      )
    },
    error = function(e) {

      # Save error message if this dataset/kernel/selection-strength fails.
      err_msg <- paste0(
        Sys.time(),
        " FAILED: dataset=", dataset_name,
        " Selection_strength=", selection_strength,
        " kernel=", kernel_mode,
        " | ", conditionMessage(e)
      )

      write(
        err_msg,
        file = file.path(
          out_dir,
          paste0(
            "results_svm_tuned_grid/svm_",
            evol_mode_comparision,
            "_grid_failures.log"
          )
        ),
        append = TRUE
      )

      return(NULL)
    }
  )

  message(Sys.time(), " Finished tune_grid()")


  ############################################################
  ## Return failed-result placeholder if tuning failed
  ############################################################
  if (is.null(svm_tuned)) {
    return(
      tibble(
        kernel = kernel_mode,
        Selection_strength = selection_strength,
        .metric = "roc_auc",
        mean = NA_real_,
        n = NA_integer_,
        std_err = NA_real_,
        .config = NA_character_,
        status = "FAILED",
        dataset = dataset_name,
        tuning_metric = "roc_auc"
      )
    )
  }


  ############################################################
  ## Save full tuning object
  ############################################################
  saveRDS(
    svm_tuned,
    file.path(
      out_dir,
      paste0(
        "results_svm_tuned_grid/svm_tuned_",
        evol_mode_comparision, "_",
        kernel_mode, "_",
        selection_strength, "_",
        dataset_name,
        ".rds"
      )
    )
  )


  ############################################################
  ## Save tuning notes
  ############################################################
  # show_notes() is useful for diagnosing failed parameter
  # combinations during tuning.
  notes_file <- file.path(
    out_dir,
    paste0(
      "results_svm_tuned_grid/svm_notes_",
      evol_mode_comparision, "_",
      kernel_mode, "_",
      selection_strength, "_",
      dataset_name,
      ".txt"
    )
  )

  writeLines(
    capture.output(show_notes(svm_tuned)),
    con = notes_file
  )


  ############################################################
  ## Select and save best hyperparameters
  ############################################################
  # Although tune_grid() calculates several metrics, this line defines
  # the metric used for model selection.
  #
  # Here, the best hyperparameter combination is chosen by maximizing
  # ROC AUC.
  best_params <- select_best(
    svm_tuned,
    metric = "roc_auc"
  )

  write.table(
    best_params,
    file = file.path(
      out_dir,
      paste0(
        "results_svm_tuned_grid/svm_best_params_",
        evol_mode_comparision, "_",
        kernel_mode, "_",
        selection_strength, "_",
        dataset_name,
        ".txt"
      )
    ),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )


  ############################################################
  ## Collect and return metrics
  ############################################################
  # This returns the full metrics table for all parameter combinations,
  # not only the selected best combination.
  results <- svm_tuned %>%
    collect_metrics() %>%
    mutate(
      kernel = kernel_mode,
      Selection_strength = selection_strength,
      dataset = dataset_name,
      status = "OK",
      tuning_metric = "roc_auc"
    )

  message(
    Sys.time(),
    " Finished SVM run: Selection_strength = ", selection_strength,
    " | kernel = ", kernel_mode,
    " | dataset = ", dataset_name
  )

  return(results)
}