### Dependencies ###
import numpy as np
import pandas as pd
import sys
### sys.argv[1] = the admixture proportion ###
### sys.argv[2] = sample size ###
### sys.argv[3] = replicate id ###
### sys.argv[4] = introgressed chromosome threshold ###


# Intialize simulation parameters.
f = float(sys.argv[1])
n = int(sys.argv[2])
rep = int(sys.argv[3])
chrom_threshold = float(sys.argv[4])

# Load the introgressed tracts.
tract_df = pd.read_csv(f'./round3_reviews/tracts/iua_f{f}_n{n}_rep{rep}.csv.gz')
# Load the windowed D and D+ results.
wind_df = pd.read_csv(f'./round3_reviews/intro/iua_f{f}_n{n}_rep{rep}.csv.gz')

# Intialize the minimum number of chromosomes.
min_chroms = int(n * 2 * chrom_threshold)
# Filter for introgressed tracts in the minimun number of chromosomes.
intro_df = tract_df[tract_df['n_intro'] >= min_chroms]
# Intialize the left window positions.
wind_lefts = wind_df.left.values
# Intialize the left and right introgressed tract positions.
intro_lefts = intro_df.left.values
intro_rights = intro_df.right.values
# Intialize the overlap.
overlaps = np.zeros(wind_lefts.size)

# For every left window position.
for i, wind_left in enumerate(wind_lefts):
    # Intialize the right position.
    wind_right = wind_left + 50_000
    # Intialize the overlap.
    overlap = 0
    # For every introgressed tract.
    for j in range(intro_lefts.size):
        # Update the amount of overlap.
        overlap += max(0, min(wind_right, intro_rights[j]) - max(wind_left, intro_lefts[j]))
    # Update the overlap.
    overlaps[i] = overlap

# Update the dataframe.
wind_df['overlap_bp'] = overlaps.astype('int')
wind_df['overlap_pc'] = overlaps / 50_000

# Export the window information.
wind_df.to_csv(f'./round3_reviews/overlap/{n}/{chrom_threshold}/iua_f{f}_n{n}_rep{rep}.csv.gz', index=False)