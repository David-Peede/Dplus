### Dependencies ###
import numpy as np
import sys
### sys.argv[1] = replicate ID ###
### sys.argv[2] = error population ###


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

# Define bootstrapping function.
def genotype_matrix_bootstrap(
    genotype_matrix,
    variable_positions,
    n_replicates,
):
    """
    ###########################################################################
    INPUT
        genotype_matrix: Simulated genome.
        variable_positions: Variable positions associated with the simulated
                            genome.
        n_replicates: Number of bootstrap replicates.
    ---------------------------------------------------------------------------
    OUTPUT: List of site pattern counts from bootstrapped replicates.
    ###########################################################################
    """
    # Intialize arrays to store bootstrapped values.
    abba_array = np.array([])
    baba_array = np.array([])
    bbaa_array = np.array([])
    baaa_array = np.array([])
    abaa_array = np.array([])
    aaba_array = np.array([])
    # For each bootstrap replicate.
    for i in range(n_replicates):
        # Intialize an empty array to store the bootstrapped genome.
        bootstrapped_genome = np.empty((0, 3))
        # For each bootstrap window (1000 x 100 kb).
        for j in range(1000):
            # Randomly generate a start and end position.
            start = np.random.randint(99_900_001)
            end = start + 100_000
            # Identify the variants that fall within the 100 kb window.
            variants = np.where(((start <= variable_positions) & (variable_positions <= end)))[0]
            # Subset the 100 kb window from the genome-wide genotype matrix.
            windowed_matrix = genotype_matrix[variants, :]
            # Concatenate 1000 windows of 100 kb to make a bootstrapped genome.
            bootstrapped_genome = np.concatenate((bootstrapped_genome, windowed_matrix))
        # Calculate site pattern counts for the bootstrapped replicate and store the results.
        abba, baba, bbaa, baaa, abaa, aaba = site_patterns(
            genotype_matrix=bootstrapped_genome,
            p1_idx=[0],
            p2_idx=[1],
            p3_idx=[2],
        )
        abba_array = np.append(abba_array, abba)
        baba_array = np.append(baba_array, baba)
        bbaa_array = np.append(bbaa_array, bbaa)
        baaa_array = np.append(baaa_array, baaa)
        abaa_array = np.append(abaa_array, abaa)
        aaba_array = np.append(aaba_array, aaba)
    return abba_array, baba_array, bbaa_array, baaa_array, abaa_array, aaba_array


# Parse comand-line arguments.
rep_id  = int(sys.argv[1])
err_pop = str(sys.argv[2])

# Load the genotype matrix and variable positions.
simulated_genotype_matrix = np.loadtxt(
    './seq_errors/geno_mats/{0}/rep_id_{1}_geno_mat.csv.gz'.format(err_pop, rep_id),
    dtype=int, delimiter=',',
)
simulated_variable_positions = np.loadtxt(
    './seq_errors/var_pos/{0}/rep_id_{1}_var_pos.csv.gz'.format(err_pop, rep_id),
    dtype=int, delimiter=',',
)

# Perform bootstrapping.
bs_abba, bs_baba, bs_bbaa, bs_baaa, bs_abaa, bs_aaba, = genotype_matrix_bootstrap(simulated_genotype_matrix, simulated_variable_positions, 1000)

# Save each bootstrapped distribution as a gzipped csv.
results_prefix = './seq_errors/bootstraps/{0}/rep_id_{1}_'.format(err_pop, rep_id)

np.savetxt(
    results_prefix+'abba.csv.gz',
    [bs_abba],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'baba.csv.gz',
    [bs_baba],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'bbaa.csv.gz',
    [bs_bbaa],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'baaa.csv.gz',
    [bs_baaa],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'abaa.csv.gz',
    [bs_abaa],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'aaba.csv.gz',
    [bs_aaba],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

