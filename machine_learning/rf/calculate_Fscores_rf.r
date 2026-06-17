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
#   Rscript calculate_Fscores_rf.r \
#       <input_file> \
#       <output_file>
#
# Example:
#   Rscript calculate_Fscores_rf.r \
#       results_all_stats_rf_cv10_over_vs_neutral_allDS.txt \
#       results_all_stats_rf_cv10_over_vs_neutral_allDS_with_fscores.txt

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop(
    "Usage: Rscript calculate_Fscores_rf.r ",
    "<input_file> <output_file>"
  )
}

input_file <- args[1]
output_file <- args[2]


############################################################
## 1. Read combined RF result file
############################################################
# Input file should already contain RF results for all datasets
# for one evolutionary comparison.
#
# Expected columns include:
#   .metric
#   mean
#   std_err
#   Selection_strength
#   dataset
#
# The dataset column should contain:
#   last_sampling
#   moments
#   autocorrelation
#   moments_autocorrelation

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

rf_all <- read.table(
  input_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)


############################################################
## 2. Check required columns
############################################################

required_cols <- c(
  ".metric",
  "mean",
  "std_err",
  "Selection_strength",
  "dataset"
)

missing_cols <- setdiff(required_cols, colnames(rf_all))

if (length(missing_cols) > 0) {
  stop(
    "Input file is missing required column(s): ",
    paste(missing_cols, collapse = ", ")
  )
}


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
## 4. Check that precision and recall are present
############################################################

required_metric_cols <- c(
  "mean_precision",
  "mean_recall"
)

missing_metric_cols <- setdiff(required_metric_cols, colnames(rf_wide))

if (length(missing_metric_cols) > 0) {
  stop(
    "Cannot calculate F-scores because the following columns are missing after pivoting: ",
    paste(missing_metric_cols, collapse = ", "),
    ". Make sure the input file contains precision and recall metrics."
  )
}


############################################################
## 5. Calculate F-scores
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
## 6. Convert F-scores back to long format
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
## 7. Combine original RF metrics and F-scores
############################################################

rf_all_with_fscores <- bind_rows(
  rf_all,
  rf_f_long
)


############################################################
## 8. Save output
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
message("Input file: ", input_file)
message("Output written to: ", output_file)