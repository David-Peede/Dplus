### Dependencies ###
import msprime
import numpy as np
import sys
### sys.argv[1] = the admixture proportion ###
### sys.argv[2] = replicate ID ###
### sys.argv[3] = sample size ###

# Define IUA model of introgression.
def iua_human_model(admix_prop):
    # Intialize demographic model.
    iua_model = msprime.Demography()
    # We assume constant and equal effective population sizes for
    # all lineages.
    iua_model.add_population(name='AFR', initial_size=10_000)
    iua_model.add_population(name='EUR', initial_size=10_000)
    iua_model.add_population(name='NEA', initial_size=10_000)
    iua_model.add_population(name='AMH', initial_size=10_000)
    iua_model.add_population(name='HUM', initial_size=10_000)
    # Introgression from the Neanderthal to the Eurasian lineage
    # occuring 1,600 generations ago with a probability of f.
    iua_model.add_mass_migration(
        time=1_600, source='EUR', dest='NEA', proportion=admix_prop,
    )
    # The African and Eurasian lineages merge into the anatomically
    # modern human lineage 4,000 generations ago.
    iua_model.add_population_split(
        time=4_000, derived=['AFR', 'EUR'], ancestral='AMH',
    )
    # The anatomically modern human and Neanderthal lineages merge
    # into the ancestral human lineage 16,000 generations ago.
    iua_model.add_population_split(
        time=16_000, derived=['AMH', 'NEA'], ancestral='HUM',
    )
    return iua_model

# Define simulator function.
def run_100mb_sim(admix_prop, seed, n):
# Define a sample size dictionary.
    n_dicc = {
        1: {'AFR': 1, 'EUR': 1, 'NEA': 1},
        100: {'AFR': 100, 'EUR': 100, 'NEA': 1},
    }
    # Simulate a 100 Mb tree-sequence.
    ts = msprime.sim_ancestry(
        samples=[
            msprime.SampleSet(n_dicc[n]['AFR'], ploidy=1, population='AFR'),
            msprime.SampleSet(n_dicc[n]['EUR'], ploidy=1, population='EUR'),
            msprime.SampleSet(n_dicc[n]['NEA'], ploidy=1, population='NEA'),
            ],
        demography=iua_human_model(admix_prop),
        sequence_length=100_000_000, recombination_rate=10**-8,
        random_seed=seed,
    )
    # Overlay mutations on the tree-sequence.
    mts = msprime.sim_mutations(
        tree_sequence=ts, rate=1.5 * 10**-8,
        model='jc69', random_seed=seed,
        discrete_genome=False,
    )
    # Extract the genotype matrix.
    genotype_matrix = mts.genotype_matrix()
    # Save the genotype matrix.
    np.savetxt(
        './sim_outputs/n_{0}/{1}/geno_mats/rep_id_{2}_geno_mat.csv.gz'.format(n, admix_prop, seed),
        genotype_matrix,
        fmt='%d',
        delimiter=',',
        )
    # Extract the variable positions.
    variable_positions = mts.tables.sites.position
    # Save the variable positions.
    np.savetxt(
        './sim_outputs/n_{0}/{1}/var_pos/rep_id_{2}_var_pos.csv.gz'.format(n, admix_prop, seed),
        [variable_positions],
        fmt='%1.15f',
        delimiter=',',
        newline='\n',
    )
    # Save the tree-sequence just for good measures.
    mts.dump('./sim_outputs/n_{0}/{1}/mut_tree_seq/rep_id_{2}_mut_tree_seq.ts'.format(n, admix_prop, seed))
    return


# Parse comand-line arguments.
f         = float(sys.argv[1])
rep_id    = int(sys.argv[2])
samp_size = int(sys.argv[3])

# Run the simulation!
run_100mb_sim(admix_prop=f, seed=rep_id, n=samp_size)
