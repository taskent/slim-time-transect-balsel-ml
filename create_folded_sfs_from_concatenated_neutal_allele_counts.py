#!/usr/bin/env python3

"""
Create a folded Site Frequency Spectrum (SFS) for BalLeRMix+.

Usage:
    python3 create_folded_sfs_from_concat_neutal_freqs_for_ballermix.py \
        <input_file> <output_file>

Input:
    A file containing allele counts per SNP (one SNP per line).
    Assumes derived allele count is in column index 2.

Output:
    Folded SFS in BalLeRMix+ format:
        derived_allele_count    sample_size    frequency

Notes:
    - Sample size is assumed to be 200 chromosomes.
    - Folding converts counts > n/2 to n - count.
    - Frequencies are normalized by total number of SNPs.
"""

import sys

# -------------------------------------------------------------------
# Command-line arguments
# -------------------------------------------------------------------
input_file = sys.argv[1]
output_file = sys.argv[2]


# -------------------------------------------------------------------
# Initialize dictionary for allele count frequencies
# -------------------------------------------------------------------
# Keys: derived allele counts (1..100)
# Values: counts of SNPs with that allele count
#
# Since sample size = 200, folded counts go from 1 to 100.
der_al_count_d = {i: 0 for i in range(1, 101)}


# -------------------------------------------------------------------
# Count total number of SNPs (lines)
# -------------------------------------------------------------------
with open(input_file, 'r') as inf:
    count_lines = len(inf.readlines())


# -------------------------------------------------------------------
# Parse input and build folded SFS
# -------------------------------------------------------------------
with open(input_file, 'r') as inf:
    for line in inf:
        l = line.strip().split()

        # Derived allele count (column 3)
        der_al_count = int(l[2])

        if der_al_count <= 100:
            # Already in minor allele range
            der_al_count_d[der_al_count] += 1

        else:
            # Fold: map to minor allele count
            folded = 200 - der_al_count
            der_al_count_d[folded] += 1


# -------------------------------------------------------------------
# Write folded SFS
# -------------------------------------------------------------------
outf = open(output_file, 'w')

for key in der_al_count_d.keys():

    # Convert counts to frequencies
    freq = der_al_count_d[key] / count_lines

    # Print to stdout (for debugging/logging)
    print(key, freq)

    # Write output
    outf.write(f'{key}\t200\t{freq}\n')


outf.close()