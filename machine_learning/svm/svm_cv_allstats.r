#!/usr/bin/env Rscript

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
## Command-line arguments
############################################################
# Usage:
#   Rscript svm_cv_allstats.r \
#       <input_file> \
#       <output_dir> \
#       <dataset_name> \
#       <evol_mode_comparison> \
#       <positive_class> \
#       <negative_class> \
#       <selection_strength> \
#       <kernel_mode> \
#       <metrics_output_file>
#
# dataset_name should be one of:
#   last_sampling
#   moments
#   autocorrelation
#   moments_autocorrelation
#
# kernel_mode should be one of:
#   linear
#   polynomial
#   radial
#
# Example:
#   Rscript svm_cv_allstats.r \
#       moments_over_vs_neutral.txt \
#       results_svm \
#       moments \
#       over_vs_neutral \
#       Overdominance \
#       Neutral \
#       2Nes=500 \
#       radial \
#       results_all_stats_svm_grid_cv_over_vs_neutral_moments.txt

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 9) {
  stop(
    "Usage: Rscript svm_cv_allstats.r ",
    "<input_file> <output_dir> <dataset_name> <evol_mode_comparison> ",
    "<positive_class> <negative_class> <selection_strength> ",
    "<kernel_mode> <metrics_output_file>"
  )
}

input_file <- args[1]
out_dir <- args[2]
dataset_name <- args[3]
evol_mode_comparison <- args[4]
positive_class <- args[5]
negative_class <- args[6]
selection_strength <- args[7]
kernel_mode <- args[8]
metrics_output_file <- args[9]


############################################################
## Create output directories
############################################################
# Main output directory.
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Directory for full tuning objects, tuning notes, best parameters,
# and failure logs.
svm_grid_dir <- file.path(out_dir, "results_svm_tuned_grid")
dir.create(svm_grid_dir, showWarnings = FALSE, recursive = TRUE)


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
## Select predictor columns
############################################################
# All columns except metadata columns are used as predictors.
#
# Required metadata columns:
#   Selection_mode
#   Selection_strength
#   Run
#
# Predictor columns:
#   all remaining summary-statistic columns

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
## Function: SVM grid tuning with grouped 10-fold CV
############################################################
# This function runs SVM classification for:
#   - one dataset type
#   - one selection-strength value
#   - one kernel type
#   - all selected summary statistics together
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
#   - controlled by positive_class and negative_class arguments
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
                            evol_mode_comparison,
                            dataset_name,
                            stat_l,
                            selection_strength,
                            kernel_mode,
                            positive_class,
                            negative_class,
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
          "results_svm_tuned_grid",
          paste0(
            "svm_",
            evol_mode_comparison,
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
      "results_svm_tuned_grid",
      paste0(
        "svm_tuned_",
        evol_mode_comparison, "_",
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
    "results_svm_tuned_grid",
    paste0(
      "svm_notes_",
      evol_mode_comparison, "_",
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
      "results_svm_tuned_grid",
      paste0(
        "svm_best_params_",
        evol_mode_comparison, "_",
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


############################################################
## Run analysis
############################################################

results <- run_svm_grid_cv(
  df = df,
  evol_mode_comparison = evol_mode_comparison,
  dataset_name = dataset_name,
  stat_l = stat_l,
  selection_strength = selection_strength,
  kernel_mode = kernel_mode,
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