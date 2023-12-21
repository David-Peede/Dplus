### Dependencies ###
import numpy as np
import sys
### sys.argv[1] = error population ###


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

# Define Danc function.
def danc(
    baaa,
    abaa,
):
    """
    ###########################################################################
    INPUT: Genome-wide BAAA and ABAA counts.
    ---------------------------------------------------------------------------
    OUTPUT: Danc value.
    ###########################################################################
    """
    # Calculate Danc.
    danc_num = (baaa - abaa)
    danc_den = (baaa + abaa)
    if (danc_den != 0):
        danc = (danc_num / danc_den)
    else:
        danc = np.nan
    return danc

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

# Define fhom function.
def fhom(
    abba,
    baba,
    aaba,
):
    """
    ###########################################################################
    INPUT: Genome-wide ABBA, BABA, and AABA counts.
    ---------------------------------------------------------------------------
    OUTPUT: fhom value.
    ---------------------------------------------------------------------------
    NOTE: This equation only works for the case of one P3 sample.
    ###########################################################################
    """
    # Calculate fhom.
    fhom_num = (abba - baba)
    fhom_den = (abba + aaba)
    if (fhom_den != 0):
        fhom = (fhom_num / fhom_den)
    else:
        fhom = np.nan
    return fhom

# Define fanc function.
def fanc(
    baaa,
    abaa,
    bbaa,
):
    """
    ###########################################################################
    INPUT: Genome-wide BAAA, ABAA, and BBAA counts.
    ---------------------------------------------------------------------------
    OUTPUT: fanc value.
    ---------------------------------------------------------------------------
    NOTE: This equation only works for the case of one P3 sample.
    ###########################################################################
    """
    # Calculate fanc.
    fanc_num = (baaa - abaa)
    fanc_den = (baaa + bbaa)
    if (fanc_den != 0):
        fanc = (fanc_num / fanc_den)
    else:
        fanc = np.nan
    return fanc

# Define f+ function.
def fplus(
    abba,
    baba,
    aaba,
    baaa,
    abaa,
    bbaa,
):
    """
    ###########################################################################
    INPUT: Genome-wide ABBA, BABA, AABA, BAAA, ABAA, and BBAA counts.
    ---------------------------------------------------------------------------
    OUTPUT: f+ value.
    ---------------------------------------------------------------------------
    NOTE: This equation only works for the case of one P3 sample.
    ###########################################################################
    """
    # Calculate f+.
    fplus_num = ((abba - baba) + (baaa - abaa))
    fplus_den = ((abba + aaba) + (baaa + bbaa))
    if (fplus_den != 0):
        fplus = (fplus_num / fplus_den)
    else:
        fplus = np.nan
    return fplus


# Parse comand-line arguments.
err_pop = str(sys.argv[1])

# Intialize arrays to store observed values.
obs_abba = np.array([])
obs_baba = np.array([])
obs_bbaa = np.array([])
obs_baaa = np.array([])
obs_abaa = np.array([])
obs_aaba = np.array([])
obs_d = np.array([])
obs_danc = np.array([])
obs_dplus = np.array([])
obs_fhom = np.array([])
obs_fanc = np.array([])
obs_fplus = np.array([])
# Calculate introgression statistics for all replicate simulations.
for rep_id in range(1, 101):
    # Load the genotype matrix.
    rep_genotype_matrix = np.loadtxt(
        './seq_errors/geno_mats/{0}/rep_id_{1}_geno_mat.csv.gz'.format(err_pop, rep_id),
        dtype=int, delimiter=',',
    )
    # Calculate site patterns.
    rep_abba, rep_baba, rep_bbaa, rep_baaa, rep_abaa, rep_aaba = site_patterns(
        genotype_matrix=rep_genotype_matrix,
        p1_idx=[0],
        p2_idx=[1],
        p3_idx=[2],
    )
    # Calculate Patterson's D.
    rep_d = pattersons_d(
        abba=rep_abba,
        baba=rep_baba,
    )
    # Calculate Danc.
    rep_danc = danc(
        baaa=rep_baaa,
        abaa=rep_abaa,
    )
    # Calculate D+.
    rep_dplus = dplus(
        abba=rep_abba,
        baba=rep_baba,
        baaa=rep_baaa,
        abaa=rep_abaa,
    )
    # Calculate fhom.
    rep_fhom = fhom(
        abba=rep_abba,
        baba=rep_baba,
        aaba=rep_aaba,
    )
    # Calculate fanc.
    rep_fanc = fanc(
        baaa=rep_baaa,
        abaa=rep_abaa,
        bbaa=rep_bbaa,
    )
    # Calculate fplus.
    rep_fplus = fplus(
        abba=rep_abba,
        baba=rep_baba,
        aaba=rep_aaba,
        baaa=rep_baaa,
        abaa=rep_abaa,
        bbaa=rep_bbaa,
    )
    # Record the observed values for this replicate.
    obs_abba = np.append(obs_abba, rep_abba)
    obs_baba = np.append(obs_baba, rep_baba)
    obs_bbaa = np.append(obs_bbaa, rep_bbaa)
    obs_baaa = np.append(obs_baaa, rep_baaa)
    obs_abaa = np.append(obs_abaa, rep_abaa)
    obs_aaba = np.append(obs_aaba, rep_aaba)
    obs_d = np.append(obs_d, rep_d)
    obs_danc = np.append(obs_danc, rep_danc)
    obs_dplus = np.append(obs_dplus, rep_dplus)
    obs_fhom = np.append(obs_fhom, rep_fhom)
    obs_fanc = np.append(obs_fanc, rep_fanc)
    obs_fplus = np.append(obs_fplus, rep_fplus)

# Save each observed result as a gzipped csv.
results_prefix = './seq_errors/obs_vals/{0}/'.format(err_pop)

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
    results_prefix+'danc.csv.gz',
    [obs_danc],
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

np.savetxt(
    results_prefix+'fhom.csv.gz',
    [obs_fhom],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'fanc.csv.gz',
    [obs_fanc],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'fplus.csv.gz',
    [obs_fplus],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)
