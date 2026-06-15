## Overview

This repository contains simulation, feature-engineering, and machine-learning pipelines used to investigate whether temporal population-genetic summary statistics improve inference of balancing selection from time-transect data.

The project focuses on two classification problems:

### 1. Distinguishing balancing selection from neutral evolution

* Overdominance
* Neutral evolution

This analysis evaluates whether balancing selection can be reliably detected from temporal population-genetic data and whether time-transect data improves model performance.

### 2. Distinguishing alternative modes of balancing selection

* Overdominance (heterozygote advantage)
* Negative frequency-dependent selection (NFDS)

This analysis evaluates whether different balancing-selection mechanisms generate distinguishable temporal genetic signatures.

For both classification tasks, machine-learning models were trained using summary statistics derived from SLiM simulations and compared across four feature sets:

* Present-day summary statistics only
* Temporal moments
* Temporal autocorrelation statistics
* Combined moments + autocorrelation statistics

Classification performance was evaluated using:

* Random Forests
* Support Vector Machines (linear, polynomial and radial kernels)

using grouped 10-fold cross-validation across simulation replicates.

All simulations were performed using SLiM.

Balancing-selection simulations were performed under eight selection-strength categories:

* 2Nes=2
* 2Nes=5
* 2Nes=10
* 2Nes=50
* 2Nes=100
* 2Nes=200
* 2Nes=500
* 2Nes=1000

where Ne denotes effective population size and s denotes the selection coefficient.

Classification performance was evaluated separately for each selection-strength category.

---

## Repository Structure

```text
.
├── simulations/
│   ├── overdominance.slim
│   ├── ndfs.slim
│   └── neutral.slim
│
├── summary_statistics/
│   ├── calculate_summary_statistics.py
│   ├── b0maf_beta1_pipeline.sh
│   ├── convert_ms2vcf.py
│   ├── create_folded_sfs_from_concatenated_neutal_allele_counts.py
│   └── moments_autocorrelation_estimation_save_datasets.r
│
├── machine_learning/
│   ├── rf/
│   │   ├── rf_cv_allstats.r
│   │   ├── rf_cv_singlestat.r
│   │   └── calculate_Fscores_rf.r
│   │
│   └── svm/
│       ├── svm_cv_allstats.r
│       ├── svm_cv_singlestat.r
│       └── extract_best_metrics_wFscores_svm.r
│
└── README.md
```

---

## Simulation Models

Three evolutionary models are implemented in SLiM:

### Neutral

No selection acts on the focal mutation.

```text
neutral.slim
```

### Overdominance

Heterozygotes have higher fitness than either homozygote.

```text
overdominance.slim
```

### Negative Frequency-Dependent Selection (NFDS)

Fitness is negatively correlated with allele frequency. Rare alleles experience higher fitness, whereas fitness decreases progressively as allele frequency increases.

```text
ndfs.slim
```

---

## Summary Statistics

The script

```bash
calculate_summary_statistics.py
```

computes population-genetic summary statistics from ms-format SLiM output.

### Usage

```bash
python3 calculate_summary_statistics.py \
    <input_ms_file> \
    <selection_mode> \
    <selection_strength> \
    <run_id> \
    <sampling_time> \
    <output_file>
```

Arguments:

| Argument           | Description                                |
| ------------------ | ------------------------------------------ |
| input_ms_file      | ms-format SLiM output file                 |
| selection_mode     | Overdominance, NFDS or Neutral             |
| selection_strength | Selection-strength category (e.g. 2Nes500) |
| run_id             | Simulation replicate identifier            |
| sampling_time      | Sampling time relative to the present      |
| output_file        | output file including summary statistics   |

Statistics include:

### Diversity statistics

* π
* Watterson's θ
* Number of singletons

### Heterozygosity statistics

* Observed heterozygosity (mean, median, maximum)
* Observed/expected heterozygosity ratio (mean, median, maximum)

### Linkage disequilibrium

* r² (median)
* Kelly's ZnS

### Haplotype statistics

* H1
* H12
* H123
* H2/H1
* Haplotype diversity
* Number of haplotypes
* iHS (median)
* nSL (median, maximum)
* EHH (mean, median)
* Pairwise differences between haplotypes (mean, median, maximum)

### Site-frequency-spectrum statistics

* Tajima's D
* Fay & Wu's H
* Fu & Li's D*
* Fu & Li's F*
* Zeng's E
* Harpending's raggedness index

### Balancing-selection-related statistics

* NCD1

### Frequency-based statistics

* Median allele frequency
* Trimmed mean allele frequency
* Fraction of intermediate-frequency variants

---

## Balancing Selection Scan Statistics

Beta1 and B0_MAF are not estimated by `calculate_summary_statistics.py`.

Instead, they are estimated using a separate pipeline consisting of:

```text
b0maf_beta1_pipeline.sh
convert_ms2vcf.py
create_folded_sfs_from_concatenated_neutal_allele_counts.py
```

### convert_ms2vcf.py

Converts ms-format SLiM output into diploid VCF format.

#### Usage

```bash
python3 convert_ms2vcf.py \
    <input_ms_file> \
    <output_vcf_file> \
    <window_size>
```

Arguments:

| Argument        | Description                                |
| --------------- | ------------------------------------------ |
| input_ms_file   | Input ms-format simulation file            |
| output_vcf_file | Output VCF file                            |
| window_size     | Genomic window centered on the target site |

---

### create_folded_sfs_from_concatenated_neutal_allele_counts.py

Constructs a folded site-frequency spectrum (SFS) from neutral simulations for use by BalLeRMix+.

#### Usage

```bash
python3 create_folded_sfs_from_concatenated_neutal_allele_counts.py \
    <input_file> \
    <output_file>
```

Arguments:

| Argument    | Description                            |
| ----------- | -------------------------------------- |
| input_file  | Concatenated neutral allele-count file |
| output_file | Folded SFS file                        |

---

### b0maf_beta1_pipeline.sh

Pipeline for estimating:

* B0_MAF (BalLeRMix+)
* Beta1 (BetaScan)

Pipeline steps:

1. Convert ms output to VCF
2. Generate BalLeRMix+ input
3. Construct folded neutral SFS
4. Run BalLeRMix+
5. Prepare BetaScan input
6. Run BetaScan

#### Usage

```bash
bash b0maf_beta1_pipeline.sh \
    <input_ms_file> \
    <output_prefix> \
    <window_size> \
    <neutral_sfs_file>
```

Arguments:

| Argument         | Description                           |
| ---------------- | ------------------------------------- |
| input_ms_file    | ms-format simulation file             |
| output_prefix    | Prefix for all output files           |
| window_size      | Window size around target site        |
| neutral_sfs_file | Folded neutral SFS used by BalLeRMix+ |

Required external software:

* BalLeRMix+
* BetaScan
* Python 3

---

## Temporal Feature Engineering

The script

```bash
moments_autocorrelation_estimation_save_datasets.r
```

generates four datasets from the time-series summary statistics.

### Usage

```bash
Rscript moments_autocorrelation_estimation_save_datasets.r \
    <input_summary_statistics_file> \
    <output_directory>
```

Arguments:

| Argument                      | Description                                        |
| ----------------------------- | -------------------------------------------------- |
| input_summary_statistics_file | Combined summary statistics table                  |
| output_directory              | Directory where processed datasets will be written |

### 1. Last Sampling

Uses summary statistics from the most recent sampling time only.

```text
dataset = last_sampling
```

### 2. Moments

Calculates:

* Mean
* Variance
* Skewness
* Kurtosis

across all sampling times.

```text
dataset = moments
```

### 3. Autocorrelation

Calculates temporal autocorrelation:

```text
lag 1–5
```

for every summary statistic.

```text
dataset = autocorrelation
```

### 4. Moments + Autocorrelation

Combines both feature sets.

```text
dataset = moments_autocorrelation
```


---

## Machine Learning Classification

Machine-learning analyses were performed using Random Forests and Support Vector Machines. For each classification task, models were trained separately for each selection-strength category and each temporal feature set.

The four supported dataset names are:

```text
last_sampling
moments
autocorrelation
moments_autocorrelation
```

The main classification comparisons are:

```text
over_vs_neutral
over_vs_nfds
```

where `over_vs_neutral` compares Overdominance against Neutral evolution, and `over_vs_nfds` compares Overdominance against Negative Frequency-Dependent Selection.

All RF and SVM scripts use grouped 10-fold cross-validation, with simulation replicate ID (`Run`) used as the grouping variable. This prevents observations from the same simulation replicate from being split across training and testing folds.

---

## Random Forest Classification

Random Forest models are implemented using `tidymodels` with the `ranger` engine.

Two types of RF analyses are included:

1. Single-statistic models
2. All-statistics models

### RF Single-Statistic Models

The script

```bash
rf_cv_singlestat.r
```

runs Random Forest classification for one summary statistic, one selection-strength category, one dataset type, and one evolutionary comparison.

For temporal datasets, the input summary-statistic name is expanded automatically depending on the dataset type:

| Dataset                   | Predictors used                                      |
| ------------------------- | ---------------------------------------------------- |
| `last_sampling`           | one predictor, e.g. `TjD`                            |
| `moments`                 | `TjD_mean`, `TjD_var`, `TjD_skewness`, `TjD_kurtosis` |
| `autocorrelation`         | `autocorrLag1_TjD` to `autocorrLag5_TjD`             |
| `moments_autocorrelation` | all moment and autocorrelation predictors            |

#### Usage

```bash
Rscript rf_cv_singlestat.r \
    <input_file> \
    <output_dir> \
    <dataset_name> \
    <evol_mode_comparison> \
    <positive_class> \
    <negative_class> \
    <summ_stat_name> \
    <selection_strength> \
    <metrics_output_file>
```

Arguments:

| Argument               | Description                                                                     |
| ---------------------- | ------------------------------------------------------------------------------- |
| `input_file`           | Input dataset generated from the temporal feature-engineering step              |
| `output_dir`           | Directory where RF prediction files and outputs will be written                 |
| `dataset_name`         | One of `last_sampling`, `moments`, `autocorrelation`, `moments_autocorrelation` |
| `evol_mode_comparison` | Name of the evolutionary comparison, e.g. `over_vs_neutral` or `over_vs_nfds`   |
| `positive_class`       | Class treated as the positive/event class, e.g. `Overdominance`                 |
| `negative_class`       | Comparison class, e.g. `Neutral` or `NFDS`                                      |
| `summ_stat_name`       | Base summary-statistic name, e.g. `TjD`, `Median_R2`, `Ncd1`                    |
| `selection_strength`   | Selection-strength category, e.g. `2Nes=500`                                    |
| `metrics_output_file`  | Output file where cross-validation metrics are appended                         |

Example:

```bash
Rscript rf_cv_singlestat.r \
    moments_over_vs_neutral.txt \
    results_rf \
    moments \
    over_vs_neutral \
    Overdominance \
    Neutral \
    TjD \
    2Nes=500 \
    results_single_stats_rf_cv10_over_vs_neutral_moments.txt
```

Outputs:

* Fold-level predictions written under:

```text
<output_dir>/rf_predictions_ind_stats/
```

* Cross-validation metrics appended to:

```text
<metrics_output_file>
```

---

### RF All-Statistics Models

The script

```bash
rf_cv_allstats.r
```

runs Random Forest classification using all available summary-statistic predictors in the input dataset.

#### Usage

```bash
Rscript rf_cv_allstats.r \
    <input_file> \
    <output_dir> \
    <dataset_name> \
    <evol_mode_comparison> \
    <positive_class> \
    <negative_class> \
    <selection_strength> \
    <metrics_output_file>
```

Arguments:

| Argument               | Description                                                                     |
| ---------------------- | ------------------------------------------------------------------------------- |
| `input_file`           | Input dataset generated from the temporal feature-engineering step              |
| `output_dir`           | Directory where RF prediction files and outputs will be written                 |
| `dataset_name`         | One of `last_sampling`, `moments`, `autocorrelation`, `moments_autocorrelation` |
| `evol_mode_comparison` | Name of the evolutionary comparison, e.g. `over_vs_neutral` or `over_vs_nfds`   |
| `positive_class`       | Class treated as the positive/event class, e.g. `Overdominance`                 |
| `negative_class`       | Comparison class, e.g. `Neutral` or `NFDS`                                      |
| `selection_strength`   | Selection-strength category, e.g. `2Nes=500`                                    |
| `metrics_output_file`  | Output file where cross-validation metrics are appended                         |

Example:

```bash
Rscript rf_cv_allstats.r \
    moments_over_vs_neutral.txt \
    results_rf \
    moments \
    over_vs_neutral \
    Overdominance \
    Neutral \
    2Nes=500 \
    results_all_stats_rf_cv10_over_vs_neutral_moments.txt
```

Outputs:

* Fold-level predictions written under:

```text
<output_dir>/rf_predictions_all_stats/
```

* Cross-validation metrics appended to:

```text
<metrics_output_file>
```

RF metrics include:

* Accuracy
* PR AUC
* Precision
* Recall
* ROC AUC

---

### RF F-score Calculation

The script

```bash
calculate_Fscores_rf.r
```

adds F-score metrics to RF output tables.

It reads RF result files for the four dataset types, combines them, and calculates:

* F1
* F2
* F0.5

from the reported precision and recall values.

#### Usage

```bash
Rscript calculate_Fscores_rf.r \
    <input_dir> \
    <evol_mode_comparison> \
    <output_file>
```

Arguments:

| Argument               | Description                                               |
| ---------------------- | --------------------------------------------------------- |
| `input_dir`            | Directory containing RF metric files                      |
| `evol_mode_comparison` | Evolutionary comparison label used in RF result filenames |
| `output_file`          | Output file containing RF metrics plus F-scores           |

Expected input filename pattern:

```text
results_all_stats_rf_cv10_<evol_mode_comparison>_<dataset>.txt
```

Example:

```bash
Rscript calculate_Fscores_rf.r \
    results_rf \
    over_vs_neutral \
    results_all_stats_rf_cv10_over_vs_neutral_allDS_with_fscores.txt
```

---

## Support Vector Machine Classification

Support Vector Machine analyses are implemented using `tidymodels` with the `kernlab` engine.

Two types of SVM analyses are included:

1. Single-statistic linear SVM models
2. All-statistics grid-tuned SVM models

---

### SVM Single-Statistic Models

The script

```bash
svm_cv_singlestat.r
```

runs a fixed linear SVM classifier for one summary statistic, one selection-strength category, one dataset type, and one evolutionary comparison.

Class weights are calculated from the class counts within the current subset and passed to `kernlab` to reduce the effect of class imbalance.

#### Usage

```bash
Rscript svm_cv_singlestat.r \
    <input_file> \
    <output_dir> \
    <dataset_name> \
    <evol_mode_comparison> \
    <positive_class> \
    <negative_class> \
    <summ_stat_name> \
    <selection_strength> \
    <metrics_output_file>
```

Arguments:

| Argument               | Description                                                                     |
| ---------------------- | ------------------------------------------------------------------------------- |
| `input_file`           | Input dataset generated from the temporal feature-engineering step              |
| `output_dir`           | Directory where SVM prediction files and outputs will be written                |
| `dataset_name`         | One of `last_sampling`, `moments`, `autocorrelation`, `moments_autocorrelation` |
| `evol_mode_comparison` | Name of the evolutionary comparison, e.g. `over_vs_neutral` or `over_vs_nfds`   |
| `positive_class`       | Class treated as the positive/event class, e.g. `Overdominance`                 |
| `negative_class`       | Comparison class, e.g. `Neutral` or `NFDS`                                      |
| `summ_stat_name`       | Base summary-statistic name, e.g. `TjD`, `Median_R2`, `Ncd1`                    |
| `selection_strength`   | Selection-strength category, e.g. `2Nes=500`                                    |
| `metrics_output_file`  | Output file where cross-validation metrics are appended                         |

Example:

```bash
Rscript svm_cv_singlestat.r \
    moments_over_vs_neutral.txt \
    results_svm \
    moments \
    over_vs_neutral \
    Overdominance \
    Neutral \
    TjD \
    2Nes=500 \
    results_single_stats_svm_cv10_over_vs_neutral_moments.txt
```

Outputs:

* Fold-level predictions written under:

```text
<output_dir>/svm_predictions_ind_stats/
```

* Cross-validation metrics appended to:

```text
<metrics_output_file>
```

---

### SVM All-Statistics Models

The script

```bash
svm_cv_allstats.r
```

runs grid-tuned SVM classification using all available predictors in the input dataset.

The supported kernels are:

```text
linear
polynomial
radial
```

For each run, the full tuning grid is evaluated using multiple metrics, and the best hyperparameter combination is selected using ROC AUC.

#### Usage

```bash
Rscript svm_cv_allstats.r \
    <input_file> \
    <output_dir> \
    <dataset_name> \
    <evol_mode_comparison> \
    <positive_class> \
    <negative_class> \
    <selection_strength> \
    <kernel_mode> \
    <metrics_output_file>
```

Arguments:

| Argument               | Description                                                                     |
| ---------------------- | ------------------------------------------------------------------------------- |
| `input_file`           | Input dataset generated from the temporal feature-engineering step              |
| `output_dir`           | Directory where SVM tuning outputs will be written                              |
| `dataset_name`         | One of `last_sampling`, `moments`, `autocorrelation`, `moments_autocorrelation` |
| `evol_mode_comparison` | Name of the evolutionary comparison, e.g. `over_vs_neutral` or `over_vs_nfds`   |
| `positive_class`       | Class treated as the positive/event class, e.g. `Overdominance`                 |
| `negative_class`       | Comparison class, e.g. `Neutral` or `NFDS`                                      |
| `selection_strength`   | Selection-strength category, e.g. `2Nes=500`                                    |
| `kernel_mode`          | One of `linear`, `polynomial`, `radial`                                         |
| `metrics_output_file`  | Output file where full tuning metrics are appended                              |

Example:

```bash
Rscript svm_cv_allstats.r \
    moments_over_vs_neutral.txt \
    results_svm \
    moments \
    over_vs_neutral \
    Overdominance \
    Neutral \
    2Nes=500 \
    radial \
    results_all_stats_svm_grid_cv_over_vs_neutral_moments.txt
```

Outputs are written under:

```text
<output_dir>/results_svm_tuned_grid/
```

This directory contains:

* Full tuning objects:

```text
svm_tuned_<evol_mode_comparison>_<kernel_mode>_<selection_strength>_<dataset_name>.rds
```

* Tuning notes:

```text
svm_notes_<evol_mode_comparison>_<kernel_mode>_<selection_strength>_<dataset_name>.txt
```

* Best hyperparameters selected by ROC AUC:

```text
svm_best_params_<evol_mode_comparison>_<kernel_mode>_<selection_strength>_<dataset_name>.txt
```

* Failure log, if any tuning runs fail:

```text
svm_<evol_mode_comparison>_grid_failures.log
```

SVM metrics include:

* Accuracy
* PR AUC
* Precision
* Recall
* ROC AUC

---

### Extracting Best SVM Models and F-scores

The script

```bash
extract_best_metrics_wFscores_svm.r
```

extracts the best SVM hyperparameter configuration for each:

```text
dataset × Selection_strength × kernel
```

The best model is selected using mean cross-validation ROC AUC. The script then extracts predictions for the ROC-AUC-selected model and calculates:

* F1
* F2
* F0.5

These F-scores are calculated for the models selected by ROC AUC; the script does not reselect models based on F-score.

#### Usage

```bash
Rscript extract_best_metrics_wFscores_svm.r \
    <input_dir> \
    <tuning_results_file> \
    <evol_mode_comparison> \
    <positive_class> \
    <negative_class> \
    <best_metrics_output> \
    <fold_fscores_output> \
    <predictions_output_dir>
```

Arguments:

| Argument                 | Description                                                               |
| ------------------------ | ------------------------------------------------------------------------- |
| `input_dir`              | Directory containing SVM tuning result files and `.rds` objects           |
| `tuning_results_file`    | Combined SVM tuning metrics file                                          |
| `evol_mode_comparison`   | Evolutionary comparison label used in SVM filenames                       |
| `positive_class`         | Class treated as the positive/event class, e.g. `Overdominance`           |
| `negative_class`         | Comparison class, e.g. `Neutral` or `NFDS`                                |
| `best_metrics_output`    | Output file containing standard metrics plus F-scores for selected models |
| `fold_fscores_output`    | Output file containing fold-level F1, F2, and F0.5 values                 |
| `predictions_output_dir` | Directory where best-model predictions will be written                    |

Example:

```bash
Rscript extract_best_metrics_wFscores_svm.r \
    results_svm/results_svm_tuned_grid \
    results_all_stats_svm_grid_cv_over_vs_neutral_allDS.txt \
    over_vs_neutral \
    Overdominance \
    Neutral \
    results_svm_best_models_over_vs_neutral_allDS.txt \
    results_svm_best_models_over_vs_neutral_allDS_foldlevel_fscores.txt \
    svm_best_predictions
```

Outputs:

* Best-model metrics table
* Fold-level F-score table
* Predictions from the ROC-AUC-selected SVM models
