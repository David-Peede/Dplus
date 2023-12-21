### Dependencies ###
import demes
import msprime
import numpy as np
import sys
### sys.argv[1] = the admixture proportion ###
### sys.argv[2] = replicate ID ###
### sys.argv[2] = sample size ###


# Define a function to run a 100 Mb simulation of the Ragsdale and Gravel 2019 model.
def run_100mb_ragsdale_sim(admix_prop, seed, n):
    # Define a sample size dictionary.
    n_dicc = {
        1: {'YRI': 1, 'CEU': 1, 'CHB': 1, 'NEA': 1},
        100: {'YRI': 100, 'CEU': 100, 'CHB': 100, 'NEA': 1},
    }
    # Load the demes graph.
    ooa_arc_intro_graph = demes.load('./yamls/ragsdale_f{0}.yaml'.format(admix_prop))
    # Simulate a 100 Mb tree-sequence.
    ts = msprime.sim_ancestry(
        samples=[
            msprime.SampleSet(n_dicc[n]['YRI'], ploidy=1, population='YRI'),
            msprime.SampleSet(n_dicc[n]['CEU'], ploidy=1, population='CEU'),
            msprime.SampleSet(n_dicc[n]['CHB'], ploidy=1, population='CHB'),
            msprime.SampleSet(n_dicc[n]['NEA'], ploidy=1, population='Neanderthal'),
            ],
        demography=msprime.Demography.from_demes(ooa_arc_intro_graph),
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
        './non_iua_models/ragsdale_2019/n_{0}/{1}/geno_mats/rep_id_{2}_geno_mat.csv.gz'.format(n, admix_prop, seed),
        genotype_matrix,
        fmt='%d',
        delimiter=',',
        )
    # Extract the variable positions.
    variable_positions = mts.tables.sites.position
    # Save the variable positions.
    np.savetxt(
        './non_iua_models/ragsdale_2019/n_{0}/{1}/var_pos/rep_id_{2}_var_pos.csv.gz'.format(n, admix_prop, seed),
        [variable_positions],
        fmt='%1.15f',
        delimiter=',',
        newline='\n',
    )
    # Save the tree-sequence just for good measures.
    mts.dump('./non_iua_models/ragsdale_2019/n_{0}/{1}/mut_tree_seq/rep_id_{2}_mut_tree_seq.ts'.format(n, admix_prop, seed))
    return


# Parse comand-line arguments.
f         = float(sys.argv[1])
rep_id    = int(sys.argv[2])
samp_size = int(sys.argv[3])

# Run the simulation!
run_100mb_ragsdale_sim(admix_prop=f, seed=rep_id, n=samp_size)
