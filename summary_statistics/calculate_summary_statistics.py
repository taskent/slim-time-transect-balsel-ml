#!/usr/bin/env python3

"""
@authors: ulas isildak (main author - ref: https://github.com/ulasisik/balancing-selection/tree/master/BaSe), ozgur taskent (added frequency based statistics)
@e-mail: isildak.ulas [at] gmail.com, ozgur.taskent86 [at] gmail.com
"""


"""
Calculate summary statistics from one ms-format SLiM output file.
Contains required modules to process simulations

Usage:
    python calculate_summary_statistics.py \
        <INPUT_FILE> <SELECTION_MODE> <SELECTION_STRENGTH> <RUN> <SAMPLING_TIME> <OUTPUT_FILE>
"""

import sys
import numpy as np
import allel


# ---------------------------------------------------------------------
# I/O helper
# ---------------------------------------------------------------------
def read_msms(filename, nchroms, seq_len):
    """
    Read an ms/msms/SLiM ms-format file into a haplotype matrix.

    Parameters
    ----------
    filename : str
        Path to ms-format file.
    nchroms : int
        Number of sampled haploid chromosomes.
    seq_len : int or float
        Simulated sequence length in bp. Positions in ms-format are scaled by this.

    Returns
    -------
    croms : np.ndarray, shape (nchroms, n_sites)
        Haplotype matrix with chromosomes as rows and segregating sites as columns.
    pos : np.ndarray, shape (n_sites,)
        Physical positions of segregating sites.
    """
    with open(filename) as handle:
        lines = handle.readlines()

    if len(lines) == 0:
        raise ValueError(f"The file {filename} is empty")

    # In ms-format, '//' marks the start of a replicate; the positions line is two lines below it.
    pointer = lines.index('//\n') + 3

    # Read and rescale positions.
    pos = lines[pointer - 1].split()
    del pos[0]  # remove the 'positions:' token
    pos = np.array(pos, dtype=float) * seq_len

    # Number of segregating sites is the length of the first haplotype line.
    n_sites = len(list(lines[pointer])) - 1

    croms = np.empty((nchroms, n_sites), dtype=object)
    for j in range(nchroms):
        row = list(lines[pointer + j])
        del row[-1]  # remove newline
        croms[j, :] = np.array(row)

    return croms.astype(int), pos


# ---------------------------------------------------------------------
# Statistic helpers
# ---------------------------------------------------------------------
def calc_median_r2(g):
    """Calculate median Rogers-Huff LD r^2 across SNP pairs."""
    gn = g.to_n_alt(fill=-1)
    ld_r = allel.rogers_huff_r(gn)
    return np.nanmedian(ld_r ** 2)


def calc_kelly_zns(g, n_pos):
    """Calculate Kelly's ZnS statistic from pairwise LD r^2 values."""
    if n_pos < 2:
        return np.nan

    gn = g.to_n_alt(fill=-1)
    ld_r = allel.rogers_huff_r(gn)
    ld_r2 = ld_r ** 2
    return (np.nansum(ld_r2) * 2.0) / (n_pos * (n_pos - 1.0))


def calc_pi(croms):
    """Calculate average pairwise differences among haplotypes."""
    n = croms.shape[0]
    if n < 2:
        return np.nan

    pairwise_diffs = []
    for i in range(n):
        for j in range(i + 1, n):
            pairwise_diffs.append(np.sum(croms[i, :] != croms[j, :]))

    return np.sum(pairwise_diffs) / ((n * (n - 1.0)) / 2.0)


def calc_faywu_h(croms):
    """Calculate Fay and Wu's H statistic."""
    n = croms.shape[0]
    counts = croms.sum(axis=0)
    s_i = np.array([np.sum(counts == i) for i in range(1, n)])
    i = np.arange(1, n)
    n_i = n - i

    theta_p = np.sum((n_i * i * s_i * 2.0) / (n * (n - 1.0)))
    theta_h = np.sum((2.0 * s_i * (i ** 2)) / (n * (n - 1.0)))
    return theta_p - theta_h


def calc_fuli_f_star(croms):
    """Calculate Fu and Li's F* statistic."""
    n = croms.shape[0]
    s = croms.shape[1]
    if s == 0:
        return np.nan

    an = np.sum(1.0 / np.arange(1, n))
    bn = np.sum(1.0 / (np.arange(1, n) ** 2))
    an1 = an + (1.0 / n)

    vfs = (((2 * n**3 + 110 * n**2 - 255 * n + 153) / (9 * n**2 * (n - 1))) +
           ((2 * (n - 1) * an) / n**2) - ((8.0 * bn) / n)) / (an**2 + bn)

    ufs = ((n / (n + 1.0) + (n + 1.0) / (3 * (n - 1.0)) -
            4.0 / (n * (n - 1.0)) +
            ((2 * (n + 1.0)) / ((n - 1.0) ** 2)) *
            (an1 - ((2.0 * n) / (n + 1.0)))) / an) - vfs

    pi_est = calc_pi(croms)
    singletons = np.sum(np.sum(croms, axis=0) == 1)
    denom = np.sqrt(ufs * s + vfs * (s ** 2.0))

    return (pi_est - (((n - 1.0) / n) * singletons)) / denom if denom != 0 else np.nan


def calc_fuli_d_star(croms):
    """Calculate Fu and Li's D* statistic."""
    n = croms.shape[0]
    s = croms.shape[1]
    if s == 0:
        return np.nan

    an = np.sum(1.0 / np.arange(1, n))
    bn = np.sum(1.0 / (np.arange(1, n) ** 2))
    an1 = an + (1.0 / n)

    cn = (2.0 * ((n * an) - 2.0 * (n - 1.0))) / ((n - 1.0) * (n - 2.0))
    dn = (cn + ((n - 2.0) / ((n - 1.0) ** 2)) +
          (2.0 / (n - 1.0)) *
          (1.5 - (2.0 * an1 - 3.0) / (n - 2.0) - 1.0 / n))

    vds = (((n / (n - 1.0)) ** 2) * bn + (an ** 2) * dn -
           (2.0 * n * an * (an + 1.0)) / ((n - 1.0) ** 2)) / (an ** 2 + bn)
    uds = ((n / (n - 1.0)) * (an - n / (n - 1.0))) - vds

    singletons = np.sum(np.sum(croms, axis=0) == 1)
    denom = np.sqrt(uds * s + vds * (s ** 2.0))

    return (((n / (n - 1.0)) * s) - (an * singletons)) / denom if denom != 0 else np.nan


def calc_zeng_e(croms):
    """Calculate Zeng's E statistic."""
    n = croms.shape[0]
    s = croms.shape[1]
    if s == 0:
        return np.nan

    an = np.sum(1.0 / np.arange(1, n))
    bn = np.sum(1.0 / (np.arange(1, n) ** 2))

    counts = croms.sum(axis=0)
    s_i = np.array([np.sum(counts == i) for i in range(1, n)])

    theta_w = s / an
    theta_l = np.sum(s_i * np.arange(1, n)) / (n - 1.0)
    theta2 = (s * (s - 1.0)) / (an ** 2 + bn)

    var1 = (n / (2.0 * (n - 1.0)) - 1.0 / an) * theta_w
    var2 = (theta2 * (bn / (an ** 2.0)) +
            2.0 * bn * (n / (n - 1.0)) ** 2.0 -
            (2.0 * (n * bn - n + 1.0)) / ((n - 1.0) * an) -
            (3.0 * n + 1.0) / (n - 1.0))

    denom = np.sqrt(var1 + var2)
    return (theta_l - theta_w) / denom if denom != 0 else np.nan


def calc_raggedness(croms):
    """Calculate Harpending's raggedness index from pairwise mismatch counts."""
    mismatches = []
    for i in range(croms.shape[0] - 1):
        for j in range(i + 1, croms.shape[0]):
            mismatches.append(np.sum(croms[i, :] != croms[j, :]))

    mismatches = np.array(mismatches)
    if mismatches.size == 0:
        return np.nan

    n_pairs = mismatches.shape[0]
    terms = []
    for i in range(1, np.max(mismatches) + 2):
        terms.append((np.sum(mismatches == i) / n_pairs -
                      np.sum(mismatches == (i - 1)) / n_pairs) ** 2)

    return np.sum(terms)


# ---------------------------------------------------------------------
# Main statistic function
# ---------------------------------------------------------------------
def sum_stats(croms, pos, nchroms, selection_mode, selection_strength, run, sampling_time):
    """
    Calculate summary statistics in the central 10 kb window.

    Central window:
        20 kb < position < 30 kb

    Frequency statistics:
        - Median_Frequency
        - Trimmed_Mean_Frequency: mean of SNP frequencies with 0.1 < f < 0.9
        - Fraction_Intermediate_Frequency: fraction of SNPs with 0.3 < f < 0.7

    """

    # -----------------------------
    # Select central 10 kb window
    # -----------------------------
    window_mask = np.logical_and(pos > 20000, pos < 30000)
    pos1 = pos[window_mask]
    croms1 = croms[:, window_mask]
    n_pos1 = croms1.shape[1]

    if n_pos1 == 0:
        raise ValueError("No segregating sites found in the central 10 kb window.")

    # Site-frequency spectrum in the central window.
    freq1 = np.asarray(croms1.sum(axis=0) / nchroms, dtype=float)

    # -----------------------------
    # Frequency-based statistics
    # -----------------------------
    median_frequency = float(np.median(freq1))

    trimmed_freq = freq1[(freq1 > 0.1) & (freq1 < 0.9)]
    trimmed_mean_frequency = float(np.mean(trimmed_freq)) if trimmed_freq.size > 0 else np.nan

    fraction_intermediate_frequency = float(np.mean((freq1 > 0.3) & (freq1 < 0.7)))

    # -----------------------------
    # Convert haplotypes for scikit-allel
    # -----------------------------
    # scikit-allel expects variants as rows and haplotypes/samples as columns for HaplotypeArray.
    haplos = croms1.T
    h1 = allel.HaplotypeArray(haplos)
    ac1 = h1.count_alleles()
    g1 = h1.to_genotypes(ploidy=2, copy=True)

    # -----------------------------
    # Diversity / SFS statistics
    # -----------------------------
 
    TjD = allel.tajima_d(ac1)
    theta_hat_w = allel.watterson_theta(pos1, ac1)

    # -----------------------------
    # Heterozygosity statistics
    # -----------------------------
    obs_het = allel.heterozygosity_observed(g1)
    af1 = ac1.to_frequencies()
    exp_het = allel.heterozygosity_expected(af1, ploidy=2)

    mean_obs_het = np.mean(obs_het)
    median_obs_het = np.median(obs_het)
    max_obs_het = np.max(obs_het)

    # Avoid division-by-zero when expected heterozygosity is zero.
    obs_over_exp_het = np.divide(
        obs_het,
        exp_het,
        out=np.zeros_like(obs_het),
        where=exp_het != 0
    )

    mean_obs_exp = np.nanmean(obs_over_exp_het)
    median_obs_exp = np.nanmedian(obs_over_exp_het)
    max_obs_exp = np.nanmax(obs_over_exp_het)

    # -----------------------------
    # LD and haplotype statistics
    # -----------------------------
    median_r2 = calc_median_r2(g1)

    garud = allel.garud_h(h1)
    h11 = garud[0]
    h12 = garud[1]
    h123 = garud[2]
    h2_h1 = garud[3]

    n_hap = np.unique(croms1, axis=0).shape[0]
    hap_div = allel.haplotype_diversity(h1)

    ehh = allel.ehh_decay(h1)
    mean_ehh = np.mean(ehh)
    median_ehh = np.median(ehh)

    ihs = allel.ihs(h1, pos1, include_edges=True)
    median_ihs = np.nanmedian(ihs)

    nsl = allel.nsl(h1)
    max_nsl = np.nanmax(nsl)
    median_nsl = np.nanmedian(nsl)

    mpd = allel.mean_pairwise_difference(ac1)
    mean_mean_pwise_dis = np.mean(mpd)
    median_mean_pwise_dis = np.median(mpd)
    max_mean_pwise_dis = np.max(mpd)

    # -----------------------------
    # Balancing Selection Statistic
    # -----------------------------
    # NCD: distance of non-fixed SNP frequencies from target frequency 0.5.
    tf = 0.5
    nonfixed_freq = freq1[freq1 < 1]
    ncd1 = np.sqrt(np.sum((nonfixed_freq - tf) ** 2) / nonfixed_freq.shape[0]) if nonfixed_freq.size > 0 else np.nan

    # -----------------------------
    # Additional neutrality / mismatch statistics
    # -----------------------------
    kellyzn = calc_kelly_zns(g1, n_pos1)
    pi_est = calc_pi(croms1)
    Hstat = calc_faywu_h(croms1)
    Ss = np.sum(np.sum(croms1, axis=0) == 1)
    Dstar = calc_fuli_d_star(croms1)
    Fstar = calc_fuli_f_star(croms1)
    ZengE = calc_zeng_e(croms1)
    rgd = calc_raggedness(croms1)

    # -----------------------------
    # Output labels and values
    # -----------------------------
    labs = [
        'Selection_mode', 'Selection_strength', 'Run', 'Sampling_time',
        'Mean(MeanPwiseDist)', 'Median(MeanPwiseDist)', 'Max(MeanPwiseDist)',
        'Tajimas D', 'Wattersons Theta',
        'Mean(ObservedHet)', 'Median(ObservedHet)', 'Max(ObservedHet)',
        'Mean(Obs/Exp Het)', 'Median(Obs/Exp Het)', 'Max(Obs/Exp Het)',
        'Median(r2)',
        'H1', 'H12', 'H123', 'H2/H1',
        'Haplotype Diversity', '# of Haplotypes',
        'Mean(EHH)', 'Median(EHH)', 'Median(ihs)',
        'Max(nsl)', 'Median(nsl)',
        'NCD1', 'KellyZns', 'pi', 'faywuH', '#ofSingletons',
        'Dstar', 'Fstar', 'ZengE', 'Raggedness',
        'Median_Frequency', 'Trimmed_Mean_Frequency', 'Fraction_Intermediate_Frequency'
    ]

    stats = [
        selection_mode, str(selection_strength), str(run), str(sampling_time), 
        mean_mean_pwise_dis, median_mean_pwise_dis, max_mean_pwise_dis,
        TjD, theta_hat_w,
        mean_obs_het, median_obs_het, max_obs_het,
        mean_obs_exp, median_obs_exp, max_obs_exp,
        median_r2,
        h11, h12, h123, h2_h1,
        hap_div, n_hap,
        mean_ehh, median_ehh, median_ihs,
        max_nsl, median_nsl,
        ncd1, kellyzn, pi_est, Hstat, Ss,
        Dstar, Fstar, ZengE, rgd,
        median_frequency, trimmed_mean_frequency, fraction_intermediate_frequency
    ]

    return labs, stats


# ---------------------------------------------------------------------
# Command-line interface
# ---------------------------------------------------------------------
def main():
    if len(sys.argv) != 7:
        raise SystemExit(
            "Usage: python calculate_summary_statistics.py "
            "<INPUT_FILE> <SELECTION_MODE> <SELECTION_STRENGTH> <RUN> <SAMPLING_TIME> <OUTPUT_FILE>"
        )

    input_file = sys.argv[1]
    selection_mode = sys.argv[2]
    selection_strength = sys.argv[3]
    run = sys.argv[4]
    sampling_time = sys.argv[5]
    output_file = sys.argv[6]

    chroms, positions = read_msms(filename=input_file, nchroms=200, seq_len=50000)
    labels, statistics = sum_stats(chroms, positions, 200, selection_mode, selection_strength, run, sampling_time)


    with open(output_file, 'a') as outf:
        print('\t'.join(str(element) for element in statistics), file=outf)


if __name__ == "__main__":
    main()
