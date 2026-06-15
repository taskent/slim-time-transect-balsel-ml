#!/usr/bin/env bash

############################################################
# Pipeline: Calculate B0_MAF and Beta1 from SLiM ms output
############################################################
#
# Usage:
#   bash b0maf_beta1_pipeline.sh <input.ms> <out_prefix> <neutral_counts_file>
#
# Example:
#   bash b0maf_beta1_pipeline.sh \
#       p_12000.ms \
#       run1_p_12000 \
#       concatenated_neutral_allele_counts.txt
#
############################################################

# Stop immediately if:
#   - any command fails
#   - an undefined variable is used
#   - any command in a pipe fails
set -euo pipefail


############################################################
# User-supplied arguments
############################################################

# ms-format SLiM output file
INPUT_MS="$1"

# Prefix used for all output files
OUT_PREFIX="$2"

# Concatenated neutral allele-count file used to build folded SFS
NEUTRAL_COUNTS="$3"


############################################################
# Fixed analysis parameters
############################################################

# Window size around target site
WINDOW=10000

# Chromosome name used in VCF / BalLeRMix+ input
CHR=1

# Recombination rate used by parse_ballermix_input.py
REC_RATE=1e-8

# BetaScan MAF threshold
BETASCAN_M=0.15


############################################################
# Derived output filenames
############################################################

# VCF converted from ms-format file
VCF="${OUT_PREFIX}.vcf"

# Parsed input file for BalLeRMix+
BALLERMIX_INPUT="${OUT_PREFIX}_ballermix_input.txt"

# Folded neutral SFS file for BalLeRMix+
SFS="${OUT_PREFIX}_neutral_sfs_folded.txt"

# BalLeRMix+ B0_MAF output
B0_OUT="${OUT_PREFIX}_B0_MAF.txt"

# Filtered BetaScan input file
BETASCAN_INPUT="${OUT_PREFIX}_betascan_input.txt"

# BetaScan Beta1 output
BETA1_OUT="${OUT_PREFIX}_Beta1.txt"


############################################################
# Step 1: Convert ms-format SLiM output to VCF
############################################################
# Keeps variants within a 10 kb window around the target site.
python3 convert_ms2vcf.py \
    "$INPUT_MS" \
    "$VCF" \
    "$WINDOW"


############################################################
# Step 2: Prepare BalLeRMix+ input from the VCF
############################################################
# Converts the VCF into the input format required by BalLeRMix+.
python3 parse_ballermix_input.py \
    --vcf "$VCF" \
    --chr "$CHR" \
    --rec_rate "$REC_RATE" \
    -o "$BALLERMIX_INPUT"


############################################################
# Step 3: Create folded Site Frequency Spectrum (SFS)
############################################################
# Builds a folded SFS from neutral simulations.
#
# Input:
#   NEUTRAL_COUNTS:
#       concatenated file of neutral allele counts
#
# Output:
#   SFS:
#       folded SFS in BalLeRMix+ format:
#           derived_allele_count   sample_size   frequency
#
# Folding:
#   counts > n/2 are mapped to n - count.
#
# This SFS is used as the neutral expectation by BalLeRMix+.
python3 create_folded_sfs_from_concatenated_neutal_allele_counts.py \
    "$NEUTRAL_COUNTS" \
    "$SFS"


############################################################
# Step 4: Run BalLeRMix+
############################################################
# Estimates B0_MAF.
#
# Important options:
#   --spect:
#       folded neutral SFS file
#
#   --usePhysPos:
#       use physical positions
#
#   --noSub:
#       no substitution information used
#
#   --MAF:
#       use minor allele frequency
#
#   --fixWinSize:
#       use fixed window size
#
#   -w:
#       window size
python3 BalLeRMix+_v1.py \
    -i "$BALLERMIX_INPUT" \
    --spect "$SFS" \
    -o "$B0_OUT" \
    --usePhysPos \
    --noSub \
    --MAF \
    --fixWinSize \
    -w "$WINDOW"


############################################################
# Step 5: Prepare BetaScan input
############################################################
# Removes the header from the BalLeRMix+ parser output.
#
# Keeps SNPs where column 3 is >= 3.
# This removes very-low-frequency variants before BetaScan.
#
# Output format:
#   chromosome   allele_count/frequency_column   position
tail -n +2 "$BALLERMIX_INPUT" | \
awk '{if($3 >= 3) {print $1"\t"$3"\t"$4}}' \
> "$BETASCAN_INPUT"


############################################################
# Step 6: Run BetaScan
############################################################
# Estimates Beta1.
#
# Options:
#   -m:
#       MAF threshold
#
#   -w:
#       window size
python3 BetaScan.py \
    -i "$BETASCAN_INPUT" \
    -m "$BETASCAN_M" \
    -w "$WINDOW" \
    -o "$BETA1_OUT"


############################################################
# Finished
############################################################
echo "Pipeline finished successfully."
echo "Input ms file:      $INPUT_MS"
echo "Output prefix:      $OUT_PREFIX"
echo "B0_MAF output:      $B0_OUT"
echo "Beta1 output:       $BETA1_OUT"