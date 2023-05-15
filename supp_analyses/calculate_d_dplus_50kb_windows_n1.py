### Dependencies ###
import numpy as np
import sys
### sys.argv[1] = error population ###
### sys.argv[2] = replicate ID ###

# Define derived allele frequency function.
def derived_allele_freq(
    genotype_matrix,
    taxon,
):
    """
    ###########################################################################
    INPUT
        genotype_matrix: Simulated genotypes from msprime.
        taxon: Lineage to calculate the derived allele frequencies for.
    ---------------------------------------------------------------------------
    OUTPUT: Derived allele frequency array for the lineage of interest.
    ---------------------------------------------------------------------------
    NOTE: This function will return counts when the sample size is 1.
    ###########################################################################
    """
    # Calculate derived allele frequencies.
    freq_array = genotype_matrix[:, taxon].sum(axis=1)/float(len(taxon))
    return freq_array

# Define site pattern function.
def site_patterns(
    genotype_matrix,
    p1_idx=[0],
    p2_idx=[1],
    p3_idx=[2],
):
    """
    ###########################################################################
    INPUT
        genotype_matrix: Simulated genotypes from msprime.
        p1_idx: P1 lineage.
        p2_idx: P2 lineage.
        p3_idx: P3 lineage.
    ---------------------------------------------------------------------------
    OUTPUT: Genome counts of ABBA, BABA, BBAA, BAAA, ABAA, and AABA sites.
    ###########################################################################
    """
    # Calculate derived allele frequencies per taxa.
    p1 = derived_allele_freq(genotype_matrix, p1_idx)
    p2 = derived_allele_freq(genotype_matrix, p2_idx)
    p3 = derived_allele_freq(genotype_matrix, p3_idx)
    # Calculate site pattern counts.
    abba = ((1 - p1) * (p2) * (p3)).sum()
    baba = ((p1) * (1 - p2) * (p3)).sum()
    bbaa = ((p1) * (p2) * (1 - p3)).sum()
    baaa = ((p1) * (1 - p2) * (1 - p3)).sum()
    abaa = ((1 - p1) * (p2) * (1 - p3)).sum()
    aaba = ((1 - p1) * (1 - p2) * (p3)).sum()
    return abba, baba, bbaa, baaa, abaa, aaba

# Define Patterson's D function.
def pattersons_d(
    abba,
    baba,
):
    """
    ###########################################################################
    INPUT: Genome-wide ABBA and BABA counts.
    ---------------------------------------------------------------------------
    OUTPUT: Patterson's D value.
    ###########################################################################
    """
    # Calculate Patterson's D.
    d_num = (abba - baba)
    d_den = (abba + baba)
    if (d_den != 0):
        d = (d_num / d_den)
    else:
        d = np.nan
    return d


# Define D+ function.
def dplus(
    abba,
    baba,
    baaa,
    abaa,
):
    """
    ###########################################################################
    INPUT: Genome-wide ABBA, BABA, BAAA, and ABAA counts.
    ---------------------------------------------------------------------------
    OUTPUT: D+ value.
    ###########################################################################
    """
    # Calculate D+.
    dplus_num = ((abba - baba) + (baaa - abaa))
    dplus_den = ((abba + baba) + (baaa + abaa))
    if (dplus_den != 0):
        dplus = (dplus_num / dplus_den)
    else:
        dplus = np.nan
    return dplus


# Parse comand-line arguments.
rep_id  = int(sys.argv[1])

# Load the genotype matrix and positions arrays.
rep_geno_mat = np.loadtxt(
    '/users/dpeede/data/data/empirical_intro_stat_benchmarking/anc-der-intro-proj/simulations/sim_outputs/n_1/0.0/geno_mats/rep_id_{0}_geno_mat.csv.gz'.format(rep_id),
    dtype=int, delimiter=',',
)
rep_var_pos = np.loadtxt(
    '/users/dpeede/data/data/empirical_intro_stat_benchmarking/anc-der-intro-proj/simulations/sim_outputs/n_1/0.0/var_pos/rep_id_{0}_var_pos.csv.gz'.format(rep_id),
    delimiter=',',
)

# Intialize arrays to store observed values.
obs_abba = np.array([])
obs_baba = np.array([])
obs_bbaa = np.array([])
obs_baaa = np.array([])
obs_abaa = np.array([])
obs_aaba = np.array([])
obs_d = np.array([])
obs_dplus = np.array([])

# For every window...
for wind in range(0, 100_000_000, 50_000):
    # Determine if there are any varibale positions in the window.
    variants = np.where(((wind <= rep_var_pos) & (rep_var_pos < (wind+50_000))))[0]
    # If there are variants to do calculations on...
    if variants.size > 0:
        # Subset the genotype matrix.
        wind_geno_mat = rep_geno_mat[variants, :]
        # Calculate site patterns.
        wind_abba, wind_baba, wind_bbaa, wind_baaa, wind_abaa, wind_aaba = site_patterns(
            genotype_matrix=wind_geno_mat,
            p1_idx=[0],
            p2_idx=[1],
            p3_idx=[2],
        )
        # Calculate Patterson's D.
        wind_d = pattersons_d(
            abba=wind_abba,
            baba=wind_baba,
        )
        # Calculate D+.
        wind_dplus = dplus(
            abba=wind_abba,
            baba=wind_baba,
            baaa=wind_baaa,
            abaa=wind_abaa,
        )
        # Record the windowed values.
        obs_abba = np.append(obs_abba, wind_abba)
        obs_baba = np.append(obs_baba, wind_baba)
        obs_bbaa = np.append(obs_bbaa, wind_bbaa)
        obs_baaa = np.append(obs_baaa, wind_baaa)
        obs_abaa = np.append(obs_abaa, wind_abaa)
        obs_aaba = np.append(obs_aaba, wind_aaba)
        obs_d = np.append(obs_d, wind_d)
        obs_dplus = np.append(obs_dplus, wind_dplus)
    # Else...
    else:
        # Record the windowed values.
        obs_abba = np.append(obs_abba, np.nan)
        obs_baba = np.append(obs_baba, np.nan)
        obs_bbaa = np.append(obs_bbaa, np.nan)
        obs_baaa = np.append(obs_baaa, np.nan)
        obs_abaa = np.append(obs_abaa, np.nan)
        obs_aaba = np.append(obs_aaba, np.nan)
        obs_d = np.append(obs_d, np.nan)
        obs_dplus = np.append(obs_dplus, np.nan)


# Save each observed result as a gzipped csv.
results_prefix = './seq_errors/wind_vals/org/rep_id_{0}_'.format(rep_id)

np.savetxt(
    results_prefix+'abba.csv.gz',
    [obs_abba],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'baba.csv.gz',
    [obs_baba],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'bbaa.csv.gz',
    [obs_bbaa],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'baaa.csv.gz',
    [obs_baaa],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'abaa.csv.gz',
    [obs_abaa],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'aaba.csv.gz',
    [obs_aaba],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'d.csv.gz',
    [obs_d],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)


np.savetxt(
    results_prefix+'dplus.csv.gz',
    [obs_dplus],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)
