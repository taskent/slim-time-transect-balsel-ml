#!/usr/bin/env python3

"""
Convert an ms-format haplotype file to a diploid VCF file.

Usage:
    python3 convert_ms2vcf.py <input.ms> <output.vcf> <window_size>

Example:
    python3 convert_ms2vcf.py p_12000.ms p_12000.vcf 10000

Notes:
    - Assumes a sequence length of 50,000 bp.
    - Assumes the selected/target site is centered at 25,000 bp.
    - Keeps only variants within the specified window around the center.
    - Assumes haplotypes are ordered such that every two rows form one diploid individual.
    - REF/ALT nucleotides are arbitrary placeholders because ms-format stores only 0/1 states.

@author: ozgur taskent
@e-mail: ozgur.taskent86 [at] gmail.com
"""

import itertools
import random
import sys


# -------------------------------------------------------------------
# Command-line arguments
# -------------------------------------------------------------------
# input_file  : ms-format input file
# output_file : VCF output file
# window_size : size of the window around 25,000 bp to retain
input_file = sys.argv[1]
output_file = sys.argv[2]
window_size = sys.argv[3]


# -------------------------------------------------------------------
# Write VCF metadata/header lines
# -------------------------------------------------------------------
# This initializes the output VCF file and writes basic VCF metadata.
# Note: sample IDs are not written here; this version directly writes
# genotype columns in the variant lines below.
with open(output_file, 'w') as outf:
    outf.write('##fileformat=VCFv4.2\n')
    outf.write('##source=convert_ms2vcf\n')
    outf.write('##FILTER=<ID=PASS,Description="All filters passed">\n')
    outf.write('##contig=<ID=1,length=50000>')
    outf.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')


# -------------------------------------------------------------------
# Read variant positions from the ms-format file
# -------------------------------------------------------------------
# ms-format stores positions as fractions between 0 and 1.
# Here, positions are converted to physical coordinates assuming
# a 50,000 bp simulated sequence.
with open(input_file, 'r') as inf:
    for line in inf:
        if line.startswith('positions'):

            # Remove the 'positions:' token and retain fractional positions
            frac_positions = line.strip().split()[1:]

            positions = []

            for frac_position in frac_positions:
                pos = float(frac_position) * 50000
                positions.append(round(pos))


# -------------------------------------------------------------------
# Extract genotypes for variants inside the target window
# -------------------------------------------------------------------
# Dictionary structure:
#     key   = physical SNP position
#     value = list of haploid allele states across sampled chromosomes
#
# For each SNP inside the requested window, the script re-opens the
# ms file and extracts the allele at that SNP index from every haplotype.
pos_genotype_d = {}

count = 0

for i, e in enumerate(positions):

    # Keep only SNPs within the window centered at 25,000 bp
    if (e >= (25000 - (int(window_size) / 2.0))) and \
       (e <= (25000 + (int(window_size) / 2.0))):

        genotype_l = []

        # Read haplotype lines and extract allele at SNP index i
        with open(input_file, 'r') as inf:
            for line in inf:

                # Haplotype lines in ms-format start with 0 or 1
                if line.startswith('0') or line.startswith('1'):
                    genotype = line.strip()[i]
                    genotype_l.append(genotype)

                # Print non-data lines that are not recognized as positions/segregating-site lines
                elif (not line.startswith('0')) and \
                     (not line.startswith('1')) and \
                     (not line.startswith('pos')) and \
                     (not line.startswith('seg')):
                    print(line)

        # Store haploid genotypes for this position
        pos_genotype_d[e] = genotype_l
        count += 1

    else:
        continue


# -------------------------------------------------------------------
# Convert haploid allele calls to diploid VCF genotypes
# -------------------------------------------------------------------
# Every two haplotypes are paired into one diploid individual:
#     hap1, hap2 -> hap1/hap2
#
# REF and ALT nucleotides are assigned randomly as complementary bases
# because ms-format does not contain nucleotide identities.
outf = open(output_file, 'a')

nucs = ['A', 'T', 'G', 'C']

for pos in pos_genotype_d.keys():

    genotype_l = pos_genotype_d[pos]

    ind_genotype_l = []
    ind_geno_l = []

    count = 0

    # Pair consecutive haplotypes into diploid genotypes
    for genotype in genotype_l:
        count += 1

        if count % 2 == 0:
            ind_geno_l.append(genotype)

            # Create diploid genotype, e.g. 0/0, 0/1, 1/1
            ind_genotype = ind_geno_l[0] + '/' + ind_geno_l[1]
            ind_genotype_l.append(ind_genotype)

            # Reset temporary genotype holder
            ind_geno_l = []

        elif count % 2 == 1:
            ind_geno_l.append(genotype)

    # Randomly assign an ALT nucleotide
    nuc_a = random.choice(nucs)

    # Assign REF as the complementary nucleotide
    if nuc_a == 'A':
        nuc_r = 'T'
    elif nuc_a == 'T':
        nuc_r = 'A'
    elif nuc_a == 'G':
        nuc_r = 'C'
    elif nuc_a == 'C':
        nuc_r = 'G'

    # Standard VCF columns:
    # CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT
    info_columns = ['1', str(pos), '.', nuc_r, nuc_a, '.', 'PASS', '.', 'GT']

    # Write one VCF line for this SNP
    print('\t'.join(info_columns + ind_genotype_l), file=outf)


# -------------------------------------------------------------------
# Close output file
# -------------------------------------------------------------------
outf.close()
