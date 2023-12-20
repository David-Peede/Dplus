### Dependencies ###
import numpy as np
import sys
### sys.argv[1] = replicate ID ###


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
    return np.array([abba, baba, bbaa, baaa, abaa, aaba])

# Define a function to make continuos positions discrete.
def modify_positions_array(org_pos_array):
    # For ever index in the positions array...
    for pos_idx in range(org_pos_array.size):
        # Round the current position.
        rounded_pos = round(org_pos_array[pos_idx])
        # Determine the indicies of the rounded position.
        rounded_pos_idx = np.where(org_pos_array.round() == rounded_pos)[0]
        # If the rounded position has multiple occurances...
        if (rounded_pos_idx.size > 1):
            # Generate the scaling array.
            scale = np.arange(rounded_pos_idx.size)
            # For every duplicated position index...
            for dup_idx in range(rounded_pos_idx.size):
                # Scale the original positions array.
                org_pos_array[rounded_pos_idx[dup_idx]] = org_pos_array[rounded_pos_idx[dup_idx]] + scale[dup_idx]
    # Round the positions array.
    round_pos_array = org_pos_array.round()
    return round_pos_array

# Define a function to induce sequencing errors.
def sim_sequencing_errors(geno_mat, mod_pos_array, seed, error_rate=1 * 10**-4, seq_len=100_000_000):
    p1_error_mat = geno_mat
    p2_error_mat = geno_mat
    p1_p2_error_mat = geno_mat
    # Determine the number of sequencing errors.
    seq_errors = round(error_rate * seq_len)
    # Determine the sites where errors will occur.
    p1_error_sites = np.random.randint(1, seq_len, seq_errors)
    p2_error_sites = np.random.randint(1, seq_len, seq_errors)
    # Concatenate the original and error sites.
    p1_all_sites = np.union1d(p1_error_sites, mod_pos_array)
    p2_all_sites = np.union1d(p2_error_sites, mod_pos_array)
    p1_p2_all_sites = np.union1d(p1_all_sites, p2_all_sites)
    # Intialize empty matricies.
    p1_error_mat = np.empty((p1_all_sites.size, 3), dtype=int)
    p2_error_mat = np.empty((p2_all_sites.size, 3), dtype=int)
    p1_p2_error_mat = np.empty((p1_p2_all_sites.size, 3), dtype=int)
    # Define a lambda function for errors.
    seq_error = lambda site: abs(site - 1)
    # For all P1 sites...
    for idx in range(p1_all_sites.size):
        # Grab the error site.
        site = p1_all_sites[idx]
        # If the error should modify an original site...
        if (site in p1_error_sites) & (site in mod_pos_array):
            # Determine if the index is an original position.
            error_idx = np.where(site == mod_pos_array)[0]
            # Grab the original site.
            org_site = geno_mat[error_idx][0]
            # Create the error.
            org_site[0] = seq_error(org_site[0])
            # Update the genotype matrix.
            p1_error_mat[idx, :] = org_site
        # Else-if the site is an original site.
        elif site in mod_pos_array:
            # Determine if the index is an original position.
            site_idx = np.where(site == mod_pos_array)[0]
            # Update the genotype matrix.
            p1_error_mat[idx, :] = geno_mat[site_idx][0]
        # Else...
        else:
            # Update the genotype matrix.
            p1_error_mat[idx, :] = np.array([1, 0, 0], dtype=int)
    # For all P2 sites...
    for idx in range(p2_all_sites.size):
        # Grab the error site.
        site = p2_all_sites[idx]
        # If the error should modify an original site...
        if (site in p2_error_sites) & (site in mod_pos_array):
            # Determine if the index is an original position.
            error_idx = np.where(site == mod_pos_array)[0]
            # Grab the original site.
            org_site = geno_mat[error_idx][0]
            # Create the error.
            org_site[1] = seq_error(org_site[1])
            # Update the genotype matrix.
            p2_error_mat[idx, :] = org_site
        # Else-if the site is an original site.
        elif site in mod_pos_array:
            # Determine if the index is an original position.
            site_idx = np.where(site == mod_pos_array)[0]
            # Update the genotype matrix.
            p2_error_mat[idx, :] = geno_mat[site_idx][0]
        # Else...
        else:
            # Update the genotype matrix.
            p2_error_mat[idx, :] = np.array([0, 1, 0], dtype=int)
    # For all P1 and P2 sites...
    for idx in range(p1_p2_all_sites.size):
        # Grab the error site.
        site = p1_p2_all_sites[idx]
        # If the sit is in the original positions and in P1 and P2.
        if (site in p1_error_sites) & (site in p2_error_sites) & (site in mod_pos_array):
            # Determine if the index is an original position.
            error_idx = np.where(site == mod_pos_array)[0]
            # Grab the original site.
            org_site = geno_mat[error_idx][0]
            # Create the errors.
            org_site[0] = seq_error(org_site[0])
            org_site[1] = seq_error(org_site[1])
            # Update the genotype matrix.
            p1_p2_error_mat[idx, :] = org_site
        # Else-if the site in the original positions and P1.
        elif (site in p1_error_sites) & (site in mod_pos_array):
            # Determine if the index is an original position.
            error_idx = np.where(site == mod_pos_array)[0]
            # Grab the original site.
            org_site = geno_mat[error_idx][0]
            # Create the errors.
            org_site[0] = seq_error(org_site[0])
            # Update the genotype matrix.
            p1_p2_error_mat[idx, :] = org_site
        # Else-if the site in the original positions and P2.
        elif (site in p2_error_sites) & (site in mod_pos_array):
            # Determine if the index is an original position.
            error_idx = np.where(site == mod_pos_array)[0]
            # Grab the original site.
            org_site = geno_mat[error_idx][0]
            # Create the errors.
            org_site[1] = seq_error(org_site[1])
            # Update the genotype matrix.
            p1_p2_error_mat[idx, :] = org_site
        # Else-if the site is an original site.
        elif site in mod_pos_array:
            # Determine if the index is an original position.
            site_idx = np.where(site == mod_pos_array)[0]
            # Update the genotype matrix.
            p1_p2_error_mat[idx, :] = geno_mat[site_idx][0]
        # Else-if the site is only in P1 and P2.
        elif (site in p1_error_sites) & (site in p2_error_sites):
            # Update the genotype matrix.
            p1_p2_error_mat[idx, :] = np.array([1, 1, 0], dtype=int)
        # Else-if the site is only in P1.
        elif site in p1_error_sites:
            # Update the genotype matrix.
            p1_p2_error_mat[idx, :] = np.array([1, 0, 0], dtype=int)
        # Else...
        else:
            # Update the genotype matrix.
            p1_p2_error_mat[idx, :] = np.array([0, 1, 0], dtype=int)
    # Calculate site patterns on error matrices.
    p1_error_patterns = site_patterns(p1_error_mat)
    p2_error_patterns = site_patterns(p2_error_mat)
    p1_p2_error_patterns = site_patterns(p1_p2_error_mat)
    # Save the genotype matrices.
    np.savetxt(
        './seq_errors/geno_mats/p1/rep_id_{0}_geno_mat.csv.gz'.format(seed),
        p1_error_mat,
        fmt='%d',
        delimiter=',',
        )
    np.savetxt(
        './seq_errors/geno_mats/p2/rep_id_{0}_geno_mat.csv.gz'.format(seed),
        p2_error_mat,
        fmt='%d',
        delimiter=',',
        )
    np.savetxt(
        './seq_errors/geno_mats/p1_p2/rep_id_{0}_geno_mat.csv.gz'.format(seed),
        p1_p2_error_mat,
        fmt='%d',
        delimiter=',',
        )
    # Save the variable positions.
    np.savetxt(
        './seq_errors/var_pos/p1/rep_id_{0}_var_pos.csv.gz'.format(seed),
        [p1_all_sites],
        fmt='%d',
        delimiter=',',
        )
    np.savetxt(
        './seq_errors/var_pos/p2/rep_id_{0}_var_pos.csv.gz'.format(seed),
        [p2_all_sites],
        fmt='%d',
        delimiter=',',
        )
    np.savetxt(
        './seq_errors/var_pos/p1_p2/rep_id_{0}_var_pos.csv.gz'.format(seed),
        [p1_p2_all_sites],
        fmt='%d',
        delimiter=',',
        )
    return

# Parse comand-line arguments.
rep_id = int(sys.argv[1])

# Load the genotype matrix and variable positions.
simulated_genotype_matrix = np.loadtxt(
    '/users/dpeede/data/data/empirical_intro_stat_benchmarking/anc-der-intro-proj/simulations/sim_outputs/n_1/0.0/geno_mats/rep_id_{0}_geno_mat.csv.gz'.format(rep_id),
    dtype=int, delimiter=',',
)
simulated_variable_positions = np.loadtxt(
    '/users/dpeede/data/data/empirical_intro_stat_benchmarking/anc-der-intro-proj/simulations/sim_outputs/n_1/0.0/var_pos/rep_id_{0}_var_pos.csv.gz'.format(rep_id),
    delimiter=',',
)

# Round the positions array.
modified_variable_positions = modify_positions_array(simulated_variable_positions)

# Simulate sequencing errors.
sim_sequencing_errors(simulated_genotype_matrix, modified_variable_positions, rep_id, error_rate=1 * 10**-4, seq_len=100_000_000)
