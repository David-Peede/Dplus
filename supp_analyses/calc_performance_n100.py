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

# Load the null distribution summary.
null_df = pd.read_csv(f'./round3_reviews/performance/iua_null_model_summary.csv.gz')
# Load the windowed D and D+ results.
overlap_df = pd.read_csv(f'./round3_reviews/overlap/{n}/{chrom_threshold}/iua_f{f}_n{n}_rep{rep}.csv.gz')

# Subset the null dataframe.
sub_null_df = null_df[null_df['n'] == n]
# Extract the p-values.
pvals = sub_null_df['pval'].values
# Extract the critical values.
d_cv_upr = sub_null_df['d_cv_upr'].values
d_cv_lwr = sub_null_df['d_cv_lwr'].values
dplus_cv_upr = sub_null_df['dplus_cv_upr'].values
dplus_cv_lwr = sub_null_df['dplus_cv_lwr'].values
# Extract the window values.
d_vals = overlap_df['D'].values
dplus_vals = overlap_df['D+'].values
# Extract the percent overlap.
percent_overlap = overlap_df['overlap_pc'].values

# For every p-value.
for i, pval in enumerate(pvals):
    # For every threshold.
    for threshold in [0.05, 0.1, 0.25]:
        # Determine which windows are significant.
        #d_sig = (d_vals >= d_cv_upr[i]) | (d_vals <= d_cv_lwr[i])
        #dplus_sig = (dplus_vals >= d_cv_upr[i]) | (dplus_vals <= d_cv_lwr[i])
        d_sig = (d_vals >= d_cv_upr[i])
        dplus_sig = (dplus_vals >= dplus_cv_upr[i])
        # Determine which windows have are considered introgressed.
        intro_winds = percent_overlap >= threshold
        # Intialize the results.
        d_results = d_sig.astype(object)
        dplus_results = dplus_sig.astype(object)
        # Determine the true positives.
        d_results = np.where((d_sig == True) & (intro_winds == True), 'TP', d_results)
        dplus_results = np.where((dplus_sig == True) & (intro_winds == True), 'TP', dplus_results)
        # Determine the true negatives.
        d_results = np.where((d_sig == False) & (intro_winds == False), 'TN', d_results)
        dplus_results = np.where((dplus_sig == False) & (intro_winds == False), 'TN', dplus_results)
        # Determine the false positives.
        d_results = np.where((d_sig == True) & (intro_winds == False), 'FP', d_results)
        dplus_results = np.where((dplus_sig == True) & (intro_winds == False), 'FP', dplus_results)
        # Determine the false negatives
        d_results = np.where((d_sig == False) & (intro_winds == True), 'FN', d_results)
        dplus_results = np.where((dplus_sig == False) & (intro_winds == True), 'FN', dplus_results)
        # Account for the undefinded values.
        d_results = np.where(np.isnan(d_vals), np.nan, d_results)
        dplus_results = np.where(np.isnan(dplus_vals), np.nan, dplus_results)
        # Account for the dropped windows.
        d_results = np.where(((percent_overlap < threshold) & (percent_overlap != 0)), np.nan, d_results)
        dplus_results = np.where(((percent_overlap < threshold) & (percent_overlap != 0)), np.nan, dplus_results)
        # Update the overlap table.
        overlap_df[f'D-{round(pval, 2)}-{threshold}'] = d_results
        overlap_df[f'D+-{round(pval, 2)}-{threshold}'] = dplus_results

# Export the performance information.
overlap_df.to_csv(f'./round3_reviews/performance/{n}/{chrom_threshold}/iua_f{f}_n{n}_rep{rep}.csv.gz', index=False)