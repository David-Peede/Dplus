### Dependencies ###
import msprime
import numpy as np
import pandas as pd
import sys
### sys.argv[1] = sample size ###
### sys.argv[2] = replicate id ###

# Define an IUA model of introgression.
def iua_human_model(f):
    # Intialize demographic model.
    iua_model = msprime.Demography()
    # We assume constant and equal effective population sizes for
    # all lineages.
    iua_model.add_population(name='P1', initial_size=10_000)
    iua_model.add_population(name='P2', initial_size=10_000)
    iua_model.add_population(name='P3', initial_size=10_000)
    iua_model.add_population(name='P12', initial_size=10_000)
    iua_model.add_population(name='P123', initial_size=10_000)
    # Introgression from the Neanderthal to the Eurasian lineage
    # occuring 1,600 generations ago with a probability of f.
    iua_model.add_mass_migration(
        time=1_600, source='P2', dest='P3', proportion=f,
    )
    # The African and Eurasian lineages merge into the anatomically
    # modern human lineage 4,000 generations ago.
    iua_model.add_population_split(
        time=4_000, derived=['P1', 'P2'], ancestral='P12',
    )
    # The anatomically modern human and Neanderthal lineages merge
    # into the ancestral human lineage 16,000 generations ago.
    iua_model.add_population_split(
        time=16_000, derived=['P12', 'P3'], ancestral='P123',
    )
    return iua_model


# Define derived allele frequency function.
def derived_allele_freq(
    genotype_matrix,
    taxon,
):
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
    # Calculate D+.
    dplus_num = ((abba - baba) + (baaa - abaa))
    dplus_den = ((abba + baba) + (baaa + abaa))
    if (dplus_den != 0):
        dplus = (dplus_num / dplus_den)
    else:
        dplus = np.nan
    return dplus

# Intialize simulation parameters.
n = int(sys.argv[1])
rep = int(sys.argv[2])

# If we are simulating a sample size of 100.
if n == 100:
    # Simulate.
    ts = msprime.sim_ancestry(
        samples=[
            msprime.SampleSet(100, ploidy=2, population='P1'),
            msprime.SampleSet(100, ploidy=2, population='P2'),
            msprime.SampleSet(1, ploidy=2, population='P3'),
        ],
        demography=iua_human_model(f=0),
        sequence_length=20_000_000,
        recombination_rate=1e-8,
        random_seed=rep,
    )
# Else we are simulating a sample size of 1.
else:
    # Simulate.
    ts = msprime.sim_ancestry(
        samples=[
            msprime.SampleSet(1, ploidy=1, population='P1'),
            msprime.SampleSet(1, ploidy=1, population='P2'),
            msprime.SampleSet(1, ploidy=1, population='P3'),
        ],
        demography=iua_human_model(f=0),
        sequence_length=20_000_000,
        recombination_rate=1e-8,
        random_seed=rep,
    )

# Overlay mutations.
mts = msprime.sim_mutations(
    tree_sequence=ts,
    rate=1.5e-8,
    model='jc69',
    random_seed=rep,
    discrete_genome=False,
)
# Extract the genotype matrix and varible positions.
geno_mat = mts.genotype_matrix()
var_pos = mts.tables.sites.position

# Intialize a dictionary to store the results.
df_dicc = {
    'left': [],
    'right': [],
    'abba': [],
    'baba': [],
    'bbaa': [],
    'baaa': [],
    'abaa': [],
    'aaba': [],
    'D': [],
    'D+': [],
}

# For every 50kb window.
for wind in range(0, 20_000_000, 50_000):
    # Determine if there are any varibale positions in the window.
    variants = np.where(((wind <= var_pos) & (var_pos < (wind+50_000))))[0]
    # If there are variants to do calculations on...
    if variants.size > 0:
        # Subset the genotype matrix.
        wind_geno_mat = geno_mat[variants, :]
        # Calculate site patterns.
        wind_abba, wind_baba, wind_bbaa, wind_baaa, wind_abaa, wind_aaba = site_patterns(
            genotype_matrix=wind_geno_mat,
            p1_idx=mts.samples(0),
            p2_idx=mts.samples(1),
            p3_idx=mts.samples(2),
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
        df_dicc['left'].append(wind)
        df_dicc['right'].append(wind+50_000)
        df_dicc['abba'].append(wind_abba)
        df_dicc['baba'].append(wind_baba)
        df_dicc['bbaa'].append(wind_bbaa)
        df_dicc['baaa'].append(wind_baaa)
        df_dicc['abaa'].append(wind_abaa)
        df_dicc['aaba'].append(wind_aaba)
        df_dicc['D'].append(wind_d)
        df_dicc['D+'].append(wind_dplus)
    # Else...
    else:
        # Record the windowed values.
        df_dicc['left'].append(wind)
        df_dicc['right'].append(wind+50_000)
        df_dicc['abba'].append(np.nan)
        df_dicc['baba'].append(np.nan)
        df_dicc['bbaa'].append(np.nan)
        df_dicc['baaa'].append(np.nan)
        df_dicc['abaa'].append(np.nan)
        df_dicc['aaba'].append(np.nan)
        df_dicc['D'].append(np.nan)
        df_dicc['D+'].append(np.nan)

# Convert the dictionary to a dataframe.
wind_df = pd.DataFrame(data=df_dicc)

# Export the the tree-sequence.
mts.dump(f'./round3_reviews/ts/iua_f0.0_n{n}_rep{rep}.ts')

# Export the window information.
wind_df.to_csv(f'./round3_reviews/null/iua_f0.0_n{n}_rep{rep}.csv.gz', index=False)