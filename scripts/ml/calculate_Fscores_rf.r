#!/usr/bin/env Rscript

############################################################
## Required packages
############################################################
# dplyr:
#   data manipulation and wrangling
#   used for:
#     bind_rows()
#     mutate()
#     select()
#
# tidyr:
#   reshaping data between long and wide formats
#   used for:
#     pivot_wider()
#     pivot_longer()
library(dplyr)
library(tidyr)


############################################################
## Command-line arguments
############################################################
# Usage:
#   Rscript add_rf_fscores.R \
#       <input_dir> \
#       <evol_mode_comparison> \
#       <output_file>
#
# Example:
#   Rscript add_rf_fscores.R \
#       results_rf \
#       overdominance_vs_neutral \
#       results_all_stats_rf_cv10_overdominance_vs_neutral_allDS_with_fscores.txt

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop(
    "Usage: Rscript add_rf_fscores.R ",
    "<input_dir> <evol_mode_comparison> <output_file>"
  )
}

input_dir <- args[1]
evol_mode_comparison <- args[2]
output_file <- args[3]


############################################################
## 1. Define dataset types
############################################################
# The RF results are expected to exist separately for each dataset:
#
#   last_sampling:
#       summary statistics from the final sampling time only
#
#   autocorrelation:
#       autocorrelation statistics across temporal samples
#
#   moments:
#       temporal moments across sampling times
#
#   moments_autocorrelation:
#       combined temporal moments and autocorrelation statistics
#
# These names should match the file names produced by the RF pipeline.

datasets <- c(
  "last_sampling",
  "autocorrelation",
  "moments",
  "moments_autocorrelation"
)


############################################################
## 2. Read and combine RF result files
############################################################
# Expected input filename pattern:
#
#   results_all_stats_rf_cv10_<evol_mode_comparison>_<dataset>.txt
#
# Example:
#
#   results_all_stats_rf_cv10_over_vs_neutral_last_sampling.txt

rf_all <- bind_rows(
  lapply(datasets, function(ds) {

    input_file <- file.path(
      input_dir,
      paste0(
        "results_all_stats_rf_cv10_",
        evol_mode_comparison, "_",
        ds,
        ".txt"
      )
    )

    if (!file.exists(input_file)) {
      stop("Input file not found: ", input_file)
    }

    read.table(
      input_file,
      header = TRUE,
      sep = "\t",
      stringsAsFactors = FALSE
    )
  })
)


############################################################
## 3. Convert RF metrics to wide format
############################################################
# Precision and recall are originally stored as separate rows.
# To calculate F-scores, they need to be on the same row for each:
#
#   Selection_strength × dataset
#
# This assumes that RF results contain:
#   .metric
#   mean
#   std_err
#   Selection_strength
#   dataset

rf_wide <- rf_all %>%
  dplyr::select(
    .metric,
    mean,
    std_err,
    Selection_strength,
    dataset
  ) %>%
  tidyr::pivot_wider(
    names_from = .metric,
    values_from = c(mean, std_err)
  )


############################################################
## 4. Calculate F-scores
############################################################
# F-score formula:
#
#   F_beta = (1 + beta^2) * precision * recall /
#            ((beta^2 * precision) + recall)
#
# F1:
#   beta = 1
#   precision and recall weighted equally
#
# F2:
#   beta = 2
#   recall weighted more strongly
#
# F0.5:
#   beta = 0.5
#   precision weighted more strongly
#
# These F-scores are calculated from the RF precision and recall
# values already present in the input results.

rf_f <- rf_wide %>%
  mutate(
    f1 = ifelse(
      (mean_precision + mean_recall) > 0,
      2 * mean_precision * mean_recall /
        (mean_precision + mean_recall),
      NA_real_
    ),

    f2 = ifelse(
      (4 * mean_precision + mean_recall) > 0,
      5 * mean_precision * mean_recall /
        (4 * mean_precision + mean_recall),
      NA_real_
    ),

    f0_5 = ifelse(
      (0.25 * mean_precision + mean_recall) > 0,
      1.25 * mean_precision * mean_recall /
        (0.25 * mean_precision + mean_recall),
      NA_real_
    )
  )


############################################################
## 5. Convert F-scores back to long format
############################################################
# The original RF result table is in long metric format.
# This converts f1, f2 and f0_5 into the same format so they can
# be appended to the original metrics.

rf_f_long <- rf_f %>%
  dplyr::select(
    Selection_strength,
    dataset,
    f1,
    f2,
    f0_5
  ) %>%
  tidyr::pivot_longer(
    cols = c(f1, f2, f0_5),
    names_to = ".metric",
    values_to = "mean"
  ) %>%
  mutate(
    std_err = NA_real_
  )


############################################################
## 6. Combine original RF metrics and F-scores
############################################################

rf_all_with_fscores <- bind_rows(
  rf_all,
  rf_f_long
)


############################################################
## 7. Save output
############################################################

write.table(
  rf_all_with_fscores,
  output_file,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


############################################################
## Finished
############################################################

message("Finished adding RF F-scores.")
message("Output written to: ", output_file)
