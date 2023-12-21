# Import packages.
import demes

# Define a function to generate yaml files.
def generate_ragsdale_yaml(f_total):
    """
    ###########################################################################
    INPUT: Admixture proportion to simulate.
    ---------------------------------------------------------------------------
    OUTPUT: Demes yaml file with a modified version of Ragsdale et al 2019.
    ###########################################################################
    """
    # Define the yaml file.
    yaml = '''description: The three population out-of-African model popularized by Gutenkunst et
  al. (2009) and augmented by archaic contributions to both Eurasian and African populations.
  Two archaic populations split early in human history, before the African expansion,
  and contribute to Eurasian populations (putative Neanderthal branch) and to the
  African branch (a deep diverging branch within Africa). Admixture is modeled as
  symmetric migration between the archaic and modern human branches, with contribution
  ending at a given time in the past. MODIFIED TO HAVE DISCRETE PULSE OF NEANDERTHAL ADMIXTURE.
time_units: generations
doi:
- https://doi.org/10.1371/journal.pgen.1008204
demes:
- name: YRI
  epochs:
  - end_time: 10344.827586206897
    start_size: 3600.0
  - end_time: 0
    start_size: 13900.0
- name: Neanderthal
  start_time: 19275.862068965518
  ancestors:
  - YRI
  epochs:
  - end_time: 0
    start_size: 3600.0
- name: ArchaicAFR
  start_time: 17206.896551724138
  ancestors:
  - YRI
  epochs:
  - end_time: 0
    start_size: 3600.0
- name: CEU
  start_time: 2093.103448275862
  ancestors:
  - YRI
  epochs:
  - end_time: 1241.3793103448277
    start_size: 880.0
  - end_time: 0
    start_size: 2300.0
    end_size: 10855.080951853866
- name: CHB
  start_time: 1241.3793103448277
  ancestors:
  - CEU
  epochs:
  - end_time: 0
    start_size: 650.0
    end_size: 65834.77001122756
migrations:
- start_time: 4310.3448275862065
  rate: 1.98e-05
  source: YRI
  dest: ArchaicAFR
- start_time: 4310.3448275862065
  rate: 1.98e-05
  source: ArchaicAFR
  dest: YRI
- end_time: 1241.3793103448277
  rate: 0.000522
  source: YRI
  dest: CEU
- start_time: 1241.3793103448277
  rate: 2.48e-05
  source: YRI
  dest: CEU
- end_time: 1241.3793103448277
  rate: 0.000522
  source: CEU
  dest: YRI
- start_time: 1241.3793103448277
  rate: 2.48e-05
  source: CEU
  dest: YRI
- rate: 0.000113
  source: CEU
  dest: CHB
- rate: 0.000113
  source: CHB
  dest: CEU
pulses:
- sources: [Neanderthal]
  dest: CEU
  proportions: [{0}]
  time: 1667.2413793103449
- sources: [Neanderthal]
  dest: CEU
  proportions: [{0}]
  time: 943.1034482758621
- sources: [Neanderthal]
  dest: CHB
  proportions: [{0}]
  time: 943.1034482758621'''.format(f_total/2)
    # Open the yaml file.
    yaml_file = open('./yamls/ragsdale_f{0}.yaml'.format(f_total), 'w')
    # Write the yaml.
    yaml_file.write(yaml)
    # Close the yaml file.
    yaml_file.close()
    return

# Define a list of admixture proportions.
admix_props = [
    0.0, 0.01, 0.02,
    0.03, 0.04, 0.05,
    0.06, 0.07, 0.08,
    0.09, 0.1, 0.2, 
    0.3, 0.4, 0.5,
]

# For every admixture proportion...
for f in admix_props:
    # Generate a yaml file.
    generate_ragsdale_yaml(f)
