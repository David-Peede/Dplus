# Import packages.
import msprime
import tskit
import numpy as np
import sys

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

# Define a function to disentangle the eight possible coalescent histories for the IUA model.
def iua_coal_hists(ts, tp3, tgf):
    """
    ###########################################################################
    INPUT
        ts: Simulated tree sequence.
        tp3: Time (P12, P3) diverged from P123.
        tgf: Time of the introgression event.
    ---------------------------------------------------------------------------
    OUTPUT
        0: The first coalesent event is between P1 and P2, and it is
           neutral+LS.
        1: The first coalesent event is between P1 and P2, and it is
           neutral+ILS.
        2: The first coalesent event is between P1 and P2, and it is
           introgression+ILS.
        3: The first coalesent event is between P1 and P3, and it is
           neutral+ILS.
        4: The first coalesent event is between P1 and P3, and it is
           introgression+ILS.
        5: The first coalesent event is between P2 and P3, and it is
           neutral+ILS.
        6: The first coalesent event is between P2 and P3, and it is
           introgression+LS
        7: The first coalesent event is between P2 and P3, and it is
           introgression+ILS
    ###########################################################################
    """
    # Grab the first and only tree in the tree sequence.
    tree = ts.first()
    # If the first coalescent event is between P1 and P2, e.g., indicies 0
    # and 1...
    if (tree.children(3) == (0, 1)):
        # Grab the birth time for the node of the first coalescent event.
        coal_time = tree.time(3)
        # If the first coalescent event takes place before T_{P3}.
        if (coal_time < tp3):
            # Record that the first coalescent history is LS(P1, P2).
            coal_history = 0
        # Else...
        else:
            # Dump the migrations table.
            migration_table = ts.dump_tables().migrations
            # If the first migration event is NOT the introgression event...
            if (migration_table.time[0] != tgf):
                # Record that the first coalescent history is neutral+ILS(P1, P2).
                coal_history = 1
            # Else...
            else:
                # Track the P1 and P2 population samples.
                tree = ts.first(tracked_samples=[1])
                # Extract all information needed for computing introgressed
                # segments from the migrations table.
                migration_nodes = migration_table.node
                migration_sources = migration_table.source
                migration_destinations = migration_table.dest
                migration_times = migration_table.time
                migration_left = migration_table.left
                migration_right = migration_table.right
                # Build masks for the introgressing population and the
                # introgression time.
                introgression_mask = migration_destinations == 2
                within_time_mask = migration_times[introgression_mask] == tgf
                # Get introgressing nodes for given introgression time and tree coordinates
                introgressing_nodes = migration_nodes[introgression_mask][within_time_mask][0] # Index 0 only for single trees.
                introgressing_left = migration_left[introgression_mask][within_time_mask][0] # Index 0 only for single trees.
                introgressing_right = migration_right[introgression_mask][within_time_mask][0] # Index 0 only for single trees.
                # If the tree interval is covered, partially or fully, by the migrating node then it is introgressed.
                if (tree.interval[0] >= introgressing_left and tree.interval[0] < introgressing_right)\
                    or (tree.interval[1] > introgressing_left and tree.interval[1] <= introgressing_right)\
                    or (tree.interval[0] < introgressing_left and tree.interval[1] > introgressing_right)\
                    or (tree.interval[0] == introgressing_left and tree.interval[1] == introgressing_right):
                    # If the introgressing node has a leaf in the recipient population...
                    if (tree.num_tracked_samples(introgressing_nodes) > 0):
                        # Record that the first coalesent event is between P1 and P2,
                        # and it is introgression+ILS.
                        coal_history = 2
                    # Else...
                    else:
                        # Record that the first coalesent event is between P1 and P2, and
                        # it is neutral+ILS.
                        coal_history = 1
    # Else-if the first coalescent event is between P1 and P3, e.g., indicies 0
    # and 2....
    elif (tree.children(3) == (0, 2)):
        # Track the P1 population sample.
        tree = ts.first(tracked_samples=[1])
        # Dump the migrations table.
        migration_table = ts.dump_tables().migrations
        # If the first migration event is NOT the introgression event...
        if (migration_table.time[0] != tgf):
            # Record that the first coalesent event is between P1 and P3, and
            # it is neutral+ILS.
            coal_history = 3
        # Else...
        else:
            # Extract all information needed for computing introgressed
            # segments from the migrations table.
            migration_nodes = migration_table.node
            migration_sources = migration_table.source
            migration_destinations = migration_table.dest
            migration_times = migration_table.time
            migration_left = migration_table.left
            migration_right = migration_table.right
            # Build masks for the introgressing population and the
            # introgression time.
            introgression_mask = migration_destinations == 2
            within_time_mask = migration_times[introgression_mask] == tgf
            # Get introgressing nodes for given introgression time and tree coordinates
            introgressing_nodes = migration_nodes[introgression_mask][within_time_mask][0] # Index 0 only for single trees.
            introgressing_left = migration_left[introgression_mask][within_time_mask][0] # Index 0 only for single trees.
            introgressing_right = migration_right[introgression_mask][within_time_mask][0] # Index 0 only for single trees.
            # If the tree interval is covered, partially or fully, by the migrating node then it is introgressed.
            if (tree.interval[0] >= introgressing_left and tree.interval[0] < introgressing_right)\
                or (tree.interval[1] > introgressing_left and tree.interval[1] <= introgressing_right)\
                or (tree.interval[0] < introgressing_left and tree.interval[1] > introgressing_right)\
                or (tree.interval[0] == introgressing_left and tree.interval[1] == introgressing_right):
                # If the introgressing node has a leaf in the recipient population...
                if (tree.num_tracked_samples(introgressing_nodes) > 0):
                    # Record that the first coalesent event is between P1 and P3,
                    # and it is introgression+ILS.
                    coal_history = 4
                # Else...
                else:
                    # Record that the first coalesent event is between P1 and P3, and
                    # it is neutral+ILS.
                    coal_history = 3
    # Else...
    else:
        # Track the P2 population sample.
        tree = ts.first(tracked_samples=[1])
        # Dump the migrations table.
        migration_table = ts.dump_tables().migrations
        # If the first migration event is NOT the introgression event...
        if (migration_table.time[0] != tgf):
            # Record that the first coalesent event is between P2 and P3, and
            # it is neutral+ILS.
            coal_history = 5
        # Else...
        else:
            # Grab the birth time for the node of the first coalescent event.
            coal_time = tree.time(3)
            # If the first coalescent event takes place before T_{P3}.
            if (coal_time < tp3):
                # Record that the first coalesent event is between P2 and P3,
                # and it is introgression+LS.
                coal_history = 6
            # Else...
            else:
                # Extract all information needed for computing introgressed
                # segments from the migrations table.
                migration_nodes = migration_table.node
                migration_sources = migration_table.source
                migration_destinations = migration_table.dest
                migration_times = migration_table.time
                migration_left = migration_table.left
                migration_right = migration_table.right
                # Build masks for the introgressing population and the
                # introgression time.
                introgression_mask = migration_destinations == 2
                within_time_mask = migration_times[introgression_mask] == tgf
                # Get introgressing nodes for given introgression time and tree coordinates
                introgressing_nodes = migration_nodes[introgression_mask][within_time_mask][0] # Index 0 only for single trees.
                introgressing_left = migration_left[introgression_mask][within_time_mask][0] # Index 0 only for single trees.
                introgressing_right = migration_right[introgression_mask][within_time_mask][0] # Index 0 only for single trees.
                # If the tree interval is covered, partially or fully, by the migrating node then it is introgressed.
                if (tree.interval[0] >= introgressing_left and tree.interval[0] < introgressing_right)\
                    or (tree.interval[1] > introgressing_left and tree.interval[1] <= introgressing_right)\
                    or (tree.interval[0] < introgressing_left and tree.interval[1] > introgressing_right)\
                    or (tree.interval[0] == introgressing_left and tree.interval[1] == introgressing_right):
                    # If the introgressing node has a leaf in the recipient population...
                    if (tree.num_tracked_samples(introgressing_nodes) > 0):
                        # Record that the first coalesent event is between P2 and P3,
                        # and it is introgression+ILS.
                        coal_history = 7
                    # Else...
                    else:
                        # Record that the first coalesent event is between P2 and P3, and
                        # it is neutral+ILS.
                        coal_history = 5
    return coal_history

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

# Define a simulation function.
def sim_unlinked_trees(kb, f):
    # Intialize a list of introgression scenarios.
    intro_ids = [2, 4, 6, 7]
    n_ils_ids = [1, 3, 5]
    # Intialize dictionary to store results.
    results_dicc = {'N_ILS': {}, 'INTRO': {}}
    # For every key...
    for key in results_dicc.keys():
        # Fill the dictionary.
        results_dicc[key]['abba'] = []
        results_dicc[key]['baba'] = []
        results_dicc[key]['baaa'] = []
        results_dicc[key]['abaa'] = []
        results_dicc[key]['d'] = []
        results_dicc[key]['danc'] = []
        results_dicc[key]['dplus'] = []
    # Intialize a log file.
    log_file = sys.stdout
    # For 100_000 replicates...
    for rep in range(1, 100_001):
        # Perform an ancestry simulation.
        ts = msprime.sim_ancestry(
            samples=[
                msprime.SampleSet(1, ploidy=1, population='P1'),
                msprime.SampleSet(1, ploidy=1, population='P2'),
                msprime.SampleSet(1, ploidy=1, population='P3'),
            ],
            demography=iua_human_model(f),
            sequence_length=kb*1_000,
            record_migrations=True,
            random_seed=rep,
        )
        # Determine the coalescent history.
        coal_hist = iua_coal_hists(ts, 16_000, 1_600)
        # Write to the logfile.
        log_file.write(str(rep)+'\t'+str(coal_hist)+'\n')
        # If the coalescent history involves introgression...
        if coal_hist in intro_ids:
            # Overlay mutations on the tree-sequence.
            mts = msprime.sim_mutations(
                tree_sequence=ts, rate=1.5 * 10**-8,
                model='jc69', random_seed=rep,
                discrete_genome=False,
            )
            # Extract the genotype matrix.
            genotype_matrix = mts.genotype_matrix()
            # Determine the site patterns.
            abba, baba, _, baaa, abaa, _ = site_patterns(
                genotype_matrix,
                p1_idx=[0],
                p2_idx=[1],
                p3_idx=[2],
            )
            # Calculate introgression statistics.
            results_dicc['INTRO']['abba'].append(abba)
            results_dicc['INTRO']['baba'].append(baba)
            results_dicc['INTRO']['baaa'].append(baaa)
            results_dicc['INTRO']['abaa'].append(abaa)
            results_dicc['INTRO']['d'].append(pattersons_d(abba, baba))
            results_dicc['INTRO']['danc'].append(danc(baaa, abaa))
            results_dicc['INTRO']['dplus'].append(dplus(abba, baba, baaa, abaa))
            # Save the tree-sequence just for good measures.
            mts.dump('./ils_v_intro/{0}kb/intro/{1}/mut_tree_seq/rep_id_{2}_mut_tree_seq.ts'.format(kb, f, rep))
        # Else-if the coalescent history is neutral ils.
        elif coal_hist in n_ils_ids:
            # Overlay mutations on the tree-sequence.
            mts = msprime.sim_mutations(
                tree_sequence=ts, rate=1.5 * 10**-8,
                model='jc69', random_seed=rep,
                discrete_genome=False,
            )
            # Extract the genotype matrix.
            genotype_matrix = mts.genotype_matrix()
            # Determine the site patterns.
            abba, baba, _, baaa, abaa, _ = site_patterns(
                genotype_matrix,
                p1_idx=[0],
                p2_idx=[1],
                p3_idx=[2],
            )
            # Calculate introgression statistics.
            results_dicc['N_ILS']['abba'].append(abba)
            results_dicc['N_ILS']['baba'].append(baba)
            results_dicc['N_ILS']['baaa'].append(baaa)
            results_dicc['N_ILS']['abaa'].append(abaa)
            results_dicc['N_ILS']['d'].append(pattersons_d(abba, baba))
            results_dicc['N_ILS']['danc'].append(danc(baaa, abaa))
            results_dicc['N_ILS']['dplus'].append(dplus(abba, baba, baaa, abaa))
            # Save the tree-sequence just for good measures.
            mts.dump('./ils_v_intro/{0}kb/n_ils/{1}/mut_tree_seq/rep_id_{2}_mut_tree_seq.ts'.format(kb, f, rep))
        # Else...
        else:
            # Continue to the next line.
            continue
    # For every statistic...
    for stat in ['abba', 'baba', 'baaa', 'abaa', 'd', 'danc', 'dplus']:
        # Export the stat as a numpy array.
        np.savetxt(
            './ils_v_intro/{0}kb/intro/{1}/{2}.csv.gz'.format(kb, f, stat),
            [np.array(results_dicc['INTRO'][stat])],
            fmt='%1.5f',
            delimiter=',',
            newline='\n',
        )
        np.savetxt(
            './ils_v_intro/{0}kb/n_ils/{1}/{2}.csv.gz'.format(kb, f, stat),
            [np.array(results_dicc['N_ILS'][stat])],
            fmt='%1.5f',
            delimiter=',',
            newline='\n',
        )
    return

sim_unlinked_trees(kb=int(sys.argv[1]), f=float(sys.argv[2]))
