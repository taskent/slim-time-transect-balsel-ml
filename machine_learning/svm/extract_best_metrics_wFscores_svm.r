#!/usr/bin/env Rscript

############################################################
## Required packages
############################################################
# dplyr:
#   data manipulation and wrangling
#   used for:
#     filter()
#     mutate()
#     select()
#     group_by()
#     summarise()
#     arrange()
#     left_join()
#     semi_join()
#
# purrr:
#   functional programming tools
#   used for:
#     pmap()   -> iterate over dataset × Selection_strength × kernel combinations
#     walk()   -> write prediction files without returning an object
#
# tidymodels:
#   modeling framework used to read and process tuning results
#   used for:
#     collect_predictions()
#
# yardstick:
#   performance metric calculations
#   used for:
#     f_meas()
#   (F1, F2 and F0.5 are calculated from the best-performing model)
library(dplyr)
library(purrr)
library(tidymodels)
library(yardstick)

############################################################
## Command-line arguments
############################################################
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 8) {
  stop(
    "Usage: Rscript extract_best_svm_models.R ",
    "<input_dir> <tuning_results_file> <evol_mode_comparison> ",
    "<positive_class> <negative_class> ",
    "<best_metrics_output> <fold_fscores_output> <predictions_output_dir>"
  )
}

input_dir <- args[1]
tuning_results_file <- args[2]
evol_mode_comparison <- args[3]
positive_class <- args[4]
negative_class <- args[5]
best_metrics_output <- args[6]
fold_fscores_output <- args[7]
predictions_output_dir <- args[8]


############################################################
## 1. Read combined SVM tuning results
############################################################
# This file should contain the combined collect_metrics() output
# from all SVM tuning runs.

svm_all <- read.table(
  file.path(input_dir, tuning_results_file),
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
) %>%
  filter(status == "OK")


############################################################
## 2. Select best configuration by ROC AUC
############################################################
# Best hyperparameter configuration is chosen within each:
#
#   dataset × Selection_strength × kernel
#
# where:
#
#   dataset:
#       - last_sampling
#           Summary statistics from the final sampling time only
#
#       - moments
#           Temporal moments (mean, SD, skewness, kurtosis)
#           estimated across the 13 sampling times
#
#       - autocorrelation
#           Temporal autocorrelation statistics (lags 1–5)
#           estimated across the 13 sampling times
#
#       - moments_autocorrelation
#           Combined moments and autocorrelation statistics
#
#   Selection_strength:
#       - 2Nes=2
#       - 2Nes=5
#       - 2Nes=10
#       - 2Nes=50
#       - 2Nes=100
#       - 2Nes=200
#       - 2Nes=500
#       - 2Nes=1000
#
#   kernel:
#       - linear
#           Linear support vector machine
#
#       - polynomial
#           Polynomial-kernel support vector machine
#
#       - radial
#           Radial basis function (RBF) support vector machine
#
# Selection criterion:
#   1. Highest mean ROC AUC across CV folds
#   2. If tied, smallest standard error
#
# This yields one best-performing hyperparameter combination
# per dataset × Selection_strength × kernel.


best_roc <- svm_all %>%
  filter(.metric == "roc_auc") %>%
  group_by(dataset, Selection_strength, kernel) %>%
  arrange(desc(mean), std_err) %>%
  slice(1) %>%
  ungroup()


############################################################
## 3. Keep already-computed metrics for best configurations
############################################################

standard_best_metrics <- svm_all %>%
  semi_join(
    best_roc %>%
      select(dataset, Selection_strength, kernel, .config),
    by = c("dataset", "Selection_strength", "kernel", ".config")
  )


############################################################
## 4. Extract predictions and calculate F-scores
############################################################
# F-scores are computed here because they were not included in
# the original SVM tuning metrics.
#
# F1:
#   balances precision and recall equally.
#
# F2:
#   gives more weight to recall.
#
# F0.5:
#   gives more weight to precision.

extract_best_fscores_from_rds <- function(dataset_name,
                                          selection_strength,
                                          kernel_mode,
                                          best_config) {

  input_dir_rds <- file.path(input_dir, "results_svm_tuned_grid")

  rds_file <- file.path(
    input_dir_rds,
    paste0(
      "svm_tuned_",
      evol_mode_comparison, "_",
      kernel_mode, "_",
      selection_strength, "_",
      dataset_name,
      ".rds"
    )
  )

  if (!file.exists(rds_file)) {
    warning("RDS file not found: ", rds_file)
    return(list(fscore_metrics = NULL, preds = NULL, fold_fscores = NULL))
  }

  obj <- readRDS(rds_file)

  preds <- collect_predictions(obj) %>%
    filter(.config == best_config)

  if (nrow(preds) == 0) {
    warning("No predictions found for config ", best_config, " in ", rds_file)
    return(list(fscore_metrics = NULL, preds = NULL, fold_fscores = NULL))
  }

  preds <- preds %>%
    mutate(
      Selection_mode = factor(
        Selection_mode,
        levels = c(positive_class, negative_class)
      ),
      .pred_class = ifelse(
        .data[[paste0(".pred_", positive_class)]] >= 0.5,
        positive_class,
        negative_class
      ),
      .pred_class = factor(
        .pred_class,
        levels = c(positive_class, negative_class)
      ),
      dataset = dataset_name,
      Selection_strength = selection_strength,
      kernel = kernel_mode
    )

  ############################################################
  ## Fold-level F-scores
  ############################################################

  fold_fscores <- preds %>%
    group_by(id) %>%
    group_modify(~ {
      bind_rows(
        f_meas(
          .x,
          truth = Selection_mode,
          estimate = .pred_class,
          beta = 1,
          event_level = "first"
        ) %>%
          mutate(.metric = "f1"),

        f_meas(
          .x,
          truth = Selection_mode,
          estimate = .pred_class,
          beta = 2,
          event_level = "first"
        ) %>%
          mutate(.metric = "f2"),

        f_meas(
          .x,
          truth = Selection_mode,
          estimate = .pred_class,
          beta = 0.5,
          event_level = "first"
        ) %>%
          mutate(.metric = "f0_5")
      )
    }) %>%
    ungroup() %>%
    mutate(
      dataset = dataset_name,
      Selection_strength = selection_strength,
      kernel = kernel_mode,
      .config = best_config
    )

  ############################################################
  ## Summary F-scores across folds
  ############################################################

  fscore_metrics <- fold_fscores %>%
    group_by(.metric, .estimator) %>%
    summarise(
      mean = mean(.estimate, na.rm = TRUE),
      n = sum(!is.na(.estimate)),
      std_err = sd(.estimate, na.rm = TRUE) / sqrt(n),
      .groups = "drop"
    ) %>%
    mutate(
      dataset = dataset_name,
      Selection_strength = selection_strength,
      kernel = kernel_mode,
      .config = best_config
    )

  list(
    fscore_metrics = fscore_metrics,
    preds = preds,
    fold_fscores = fold_fscores
  )
}


############################################################
## 5. Run over all best configurations
############################################################

extracted <- pmap(
  list(
    best_roc$dataset,
    best_roc$Selection_strength,
    best_roc$kernel,
    best_roc$.config
  ),
  extract_best_fscores_from_rds
)


############################################################
## 6. Combine standard metrics and F-scores
############################################################

fscore_metrics <- bind_rows(lapply(extracted, function(x) x$fscore_metrics))

svm_best_metrics <- bind_rows(
  standard_best_metrics,
  fscore_metrics
)

write.table(
  svm_best_metrics,
  best_metrics_output,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


############################################################
## 7. Save predictions for best configurations
############################################################

dir.create(
  predictions_output_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

walk(extracted, function(x) {

  if (is.null(x$preds) || nrow(x$preds) == 0) return(NULL)

  ds <- unique(x$preds$dataset)
  selection_strength <- unique(x$preds$Selection_strength)
  ker <- unique(x$preds$kernel)

  out_file <- file.path(
    predictions_output_dir,
    paste0(
      "svm_predictions_all_stats_",
      ds, "_",
      selection_strength, "_",
      ker, "_",
      evol_mode_comparison,
      ".txt"
    )
  )

  write.table(
    x$preds,
    out_file,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
})


############################################################
## 8. Save fold-level F-scores
############################################################

svm_best_fold_fscores <- bind_rows(
  lapply(extracted, function(x) x$fold_fscores)
)

write.table(
  svm_best_fold_fscores,
  fold_fscores_output,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


############################################################
## Finished
############################################################

message("Finished extracting best SVM models and F-scores.")