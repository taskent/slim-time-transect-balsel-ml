#!/usr/bin/env Rscript

############################################################
## Required packages
############################################################
# dplyr:
#   data manipulation, filtering, grouping, summarising, across()
#
# e1071:
#   skewness() and kurtosis()
library(dplyr)
library(e1071)


############################################################
## Command-line arguments
############################################################
# Usage:
#   Rscript calculate_moments_autocorrelation.R \
#       <input_summary_statistics_file> \
#       <last_sampling_bal_sel_output> \
#       <last_sampling_over_vs_neutral_output> \
#       <moments_bal_sel_output> \
#       <moments_over_vs_neutral_output> \
#       <autocorr_bal_sel_output> \
#       <autocorr_over_vs_neutral_output> \
#       <moments_autocorr_bal_sel_output> \
#       <moments_autocorr_over_vs_neutral_output>

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 9) {
  stop(
    "Usage: Rscript calculate_moments_autocorrelation.R ",
    "<input_summary_statistics_file> ",
    "<last_sampling_bal_sel_output> ",
    "<last_sampling_over_vs_neutral_output> ",
    "<moments_bal_sel_output> ",
    "<moments_over_vs_neutral_output> ",
    "<autocorr_bal_sel_output> ",
    "<autocorr_over_vs_neutral_output>",
    "<moments_autocorr_bal_sel_output> ",
    "<moments_autocorr_over_vs_neutral_output>"
  )
}

input_file <- args[1]

last_sampling_bal_sel_output <- args[2]
last_sampling_over_vs_neutral_output <- args[3]

moments_bal_sel_output <- args[4]
moments_over_vs_neutral_output <- args[5]

autocorr_bal_sel_output <- args[6]
autocorr_over_vs_neutral_output <- args[7]

moments_autocorr_bal_sel_output <- args[8]
moments_autocorr_over_vs_neutral_output <- args[9]


############################################################
## 1. Read summary-statistics file
############################################################

df_summary_stats <- read.table(
  input_file,
  header = FALSE,
  sep = "\t"
)


############################################################
## 2. Rename columns with compact R-friendly names
############################################################

colnames(df_summary_stats) <- c(
  "Selection_mode",
  "Selection_strength",
  "Run",
  "Sampling_time",

  "Mean_MeanPwiseDist",
  "Median_MeanPwiseDist",
  "Max_MeanPwiseDist",

  "TjD",
  "Theta_hat_w",

  "Mean_Obs_Het",
  "Median_Obs_Het",
  "Max_Obs_Het",

  "Mean_Obs_Exp",
  "Median_Obs_Exp",
  "Max_Obs_Exp",

  "Median_R2",

  "H1",
  "H12",
  "H123",
  "H2_H1",

  "Hap_Div",
  "N_Hap",

  "Mean_Ehh",
  "Median_Ehh",
  "Median_Ihs",

  "Max_Nsl",
  "Median_Nsl",

  "Ncd1",
  "Kellyzn",
  "Pi_est",
  "Hstat",
  "Ss",
  "Dstar",
  "Fstar",
  "ZengE",
  "Rgd",

  "Median_Freq",
  "Trim_Mean_Freq",
  "Frac_Inter_Freq",

  "Max_Beta1",
  "Max_B0maf"
)


############################################################
## Optional: column descriptions
############################################################
# This object is useful for documentation and can be copied into
# the README or supplementary material.

column_descriptions <- c(
  Selection_mode = "Selection model, i.e. Overdominance, NFDS or Neutral",
  Selection_strength = "Selection strength category, e.g. 2Nes1000",
  Run = "Simulation replicate ID",
  Sampling_time = "Sampling time",

  Mean_MeanPwiseDist = "Mean of mean pairwise differences across SNPs",
  Median_MeanPwiseDist = "Median of mean pairwise differences across SNPs",
  Max_MeanPwiseDist = "Maximum mean pairwise difference across SNPs",

  TjD = "Tajima's D",
  Theta_hat_w = "Watterson's theta",

  Mean_Obs_Het = "Mean observed heterozygosity",
  Median_Obs_Het = "Median observed heterozygosity",
  Max_Obs_Het = "Maximum observed heterozygosity",

  Mean_Obs_Exp = "Mean observed/expected heterozygosity ratio",
  Median_Obs_Exp = "Median observed/expected heterozygosity ratio",
  Max_Obs_Exp = "Maximum observed/expected heterozygosity ratio",

  Median_R2 = "Median pairwise LD r^2",

  H1 = "Garud's H1",
  H12 = "Garud's H12",
  H123 = "Garud's H123",
  H2_H1 = "Garud's H2/H1",

  Hap_Div = "Haplotype diversity",
  N_Hap = "Number of unique haplotypes",

  Mean_Ehh = "Mean EHH decay",
  Median_Ehh = "Median EHH decay",
  Median_Ihs = "Median iHS",
  Max_Nsl = "Maximum nSL",
  Median_Nsl = "Median nSL",

  Ncd1 = "Distance of allele frequencies from 0.5",
  Kellyzn = "Kelly's ZnS",
  Pi_est = "Average pairwise differences, pi",
  Hstat = "Fay and Wu's H",
  Ss = "Number of singleton SNPs",
  Dstar = "Fu and Li's D*",
  Fstar = "Fu and Li's F*",
  ZengE = "Zeng's E",
  Rgd = "Harpending's raggedness index",

  Median_Freq = "Median allele frequency in the central window",
  Trim_Mean_Freq = "Mean allele frequency after excluding rare and nearly fixed variants",
  Frac_Inter_Freq = "Fraction of intermediate-frequency variants",

  Max_Beta1 = "Maximum BetaScan Beta1 statistic",
  Max_B0maf = "Maximum BalLeRMix B0_MAF statistic"
)


############################################################
## 3. Extract last-sampling data
############################################################
# Last sampling corresponds to Sampling_time == 0.
# Two classification datasets are saved:
#   1. balancing selection comparison:
#        Overdominance vs NFDS
#   2. balancing selection vs neutral comparison:
#        Overdominance vs Neutral

df_last_sampling_summary_stats <- df_summary_stats %>%
  filter(Sampling_time == 0)

df_last_sampling_summary_stats_bal_sel <- df_last_sampling_summary_stats %>%
  filter(Selection_mode %in% c("Overdominance", "NFDS"))

df_last_sampling_summary_stats_overVSneutral <- df_last_sampling_summary_stats %>%
  filter(Selection_mode %in% c("Overdominance", "Neutral"))


############################################################
## 4. Estimate temporal moments for each statistic
############################################################
# For each Selection_mode × Selection_strength × Run,
# calculate mean, SD, skewness, and kurtosis across sampling times.
#
# na.rm = TRUE:
#   missing values are ignored for moment estimation.

df_moments_summary_stats <- df_summary_stats %>%
  group_by(Selection_mode, Selection_strength, Run) %>%
  summarise(
    across(
      Mean_MeanPwiseDist:Max_B0maf,
      list(
        mean     = ~ mean(.x, na.rm = TRUE),
        var      = ~ var(.x, na.rm = TRUE),
        skewness = ~ skewness(.x, na.rm = TRUE),
        kurtosis = ~ kurtosis(.x, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  )

df_moments_summary_stats_bal_sel <- df_moments_summary_stats %>%
  filter(Selection_mode %in% c("Overdominance", "NFDS"))

df_moments_summary_stats_overVSneutral <- df_moments_summary_stats %>%
  filter(Selection_mode %in% c("Overdominance", "Neutral"))


############################################################
## 5. Helper function for autocorrelation at a specific lag
############################################################
# acf() returns lag 0 first, so:
#   lag 1 = acf[2, 1, 1]
#   lag 2 = acf[3, 1, 1]


get_acf_lag <- function(x, lag) {
  acf(
    x,
    type = "correlation",
    plot = FALSE,
    na.action = na.pass
  )$acf[lag + 1, 1, 1]
}


############################################################
## 6. Estimate autocorrelation statistics for lags 1–5
############################################################
# For each Selection_mode × Selection_strength × Run:
#   1. order observations by Sampling_time
#   2. calculate autocorrelation at lags 1–5
#      for each statistic across the temporal samples.

df_autocorr_summary_stats <- df_summary_stats %>%
  group_by(Selection_mode, Selection_strength, Run) %>%
  arrange(Sampling_time, .by_group = TRUE) %>%
  summarise(
    across(
      Mean_MeanPwiseDist:Max_B0maf,
      list(
        autocorrLag1 = ~ get_acf_lag(.x, 1),
        autocorrLag2 = ~ get_acf_lag(.x, 2),
        autocorrLag3 = ~ get_acf_lag(.x, 3),
        autocorrLag4 = ~ get_acf_lag(.x, 4),
        autocorrLag5 = ~ get_acf_lag(.x, 5)
      ),
      .names = "{.fn}_{.col}"
    ),
    .groups = "drop"
  )

df_autocorr_summary_stats_bal_sel <- df_autocorr_summary_stats %>%
  filter(Selection_mode %in% c("Overdominance", "NFDS"))

df_autocorr_summary_stats_overVSneutral <- df_autocorr_summary_stats %>%
  filter(Selection_mode %in% c("Overdominance", "Neutral"))


############################################################
## 7. Merge moments and autocorrelation  statistics 
############################################################
df_moments_autocorr_summary_stats <- inner_join(
  df_moments_summary_stats,
  df_autocorr_summary_stats,
  by = c("Selection_mode", "Selection_strength", "Run")
)

df_moments_autocorr_summary_stats_bal_sel <- df_moments_autocorr_summary_stats %>%
  filter(Selection_mode %in% c("Overdominance", "NFDS"))

df_moments_autocorr_summary_stats_overVSneutral <- df_moments_autocorr_summary_stats %>%
  filter(Selection_mode %in% c("Overdominance", "Neutral"))


############################################################
## 8. Save output files
############################################################

write.table(
  df_last_sampling_summary_stats_bal_sel,
  last_sampling_bal_sel_output,
  sep = "\t",
  col.names = TRUE,
  quote = FALSE,
  row.names = FALSE
)

write.table(
  df_last_sampling_summary_stats_overVSneutral,
  last_sampling_over_vs_neutral_output,
  sep = "\t",
  col.names = TRUE,
  quote = FALSE,
  row.names = FALSE
)

write.table(
  df_moments_summary_stats_bal_sel,
  moments_bal_sel_output,
  sep = "\t",
  col.names = TRUE,
  quote = FALSE,
  row.names = FALSE
)

write.table(
  df_moments_summary_stats_overVSneutral,
  moments_over_vs_neutral_output,
  sep = "\t",
  col.names = TRUE,
  quote = FALSE,
  row.names = FALSE
)

write.table(
  df_autocorr_summary_stats_bal_sel,
  autocorr_bal_sel_output,
  sep = "\t",
  col.names = TRUE,
  quote = FALSE,
  row.names = FALSE
)

write.table(
  df_autocorr_summary_stats_overVSneutral,
  autocorr_over_vs_neutral_output,
  sep = "\t",
  col.names = TRUE,
  quote = FALSE,
  row.names = FALSE
)

write.table(
  df_moments_autocorr_summary_stats_bal_sel,
  moments_autocorr_bal_sel_output,
  sep = "\t",
  col.names = TRUE,
  quote = FALSE,
  row.names = FALSE
)

write.table(
  df_moments_autocorr_summary_stats_overVSneutral,
  moments_autocorr_over_vs_neutral_output,
  sep = "\t",
  col.names = TRUE,
  quote = FALSE,
  row.names = FALSE
)


############################################################
## Finished
############################################################

message("Finished writing output files.")
