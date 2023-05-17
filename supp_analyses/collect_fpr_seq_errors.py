### Dependencies ###
import numpy as np
import sys
### sys.argv[1] = error population ###

# Define a function to extract observed introgression values.
def err_wind_vals(err_pop):
    # Define the file path for the results.
    results_path = './seq_errors/wind_vals/{0}/rep_id_'.format(err_pop)
    # Intialize arrays to store observed values.
    wind_d = np.array([])
    wind_dplus = np.array([])
    # For all replicates.
    for rep_id in range(1, 101):
        # Load the observed values.
        d = np.loadtxt(
            results_path+'{0}_d.csv.gz'.format(rep_id),
            delimiter=',',
        )
        dplus = np.loadtxt(
            results_path+'{0}_dplus.csv.gz'.format(rep_id),
            delimiter=',',
        )
        # Append the observed value arrays.
        wind_d = np.append(wind_d, d)
        wind_dplus = np.append(wind_dplus, dplus)
    # Create a dictionary where key is the introgression metric and its
    # corresponding numpy array is the value.
    wind_vals_dict = {
        'd': wind_d,
        'dplus': wind_dplus,
    }
    return wind_vals_dict

# Define a function to calculate power for introgression detection metrics.
def err_wind_power(err_pop, wind_dict, pval):
    # Calculate p-values.
    d_pvals = np.array([np.quantile(d, (1-np.round(pval/2, 4))) for d in wind_dict['d']])
    dplus_pvals = np.array([np.quantile(dplus, (1-np.round(pval/2, 4))) for dplus in wind_dict['dplus']])
    # Calculate the power ie the number of significant replicates out of 100.
    d_power = np.nanmean(d_pvals < pval)
    dplus_power = np.nanmean(dplus_pvals < pval)
    # Create a dictionary where key is the introgression metric and its
    # corresponding numpy array is the value.
    power_dict = {
        'd': d_power, 'dplus': dplus_power,
    }
    return power_dict

# Load the results.
winds = err_wind_vals(str(sys.argv[1]))


# Determine the critical values.
cvs = np.arange(0.01, 1.01, 0.01)
# Intialize dictionaries.
power = {'d': [], 'dplus': []}
# For all critical values...
for cv in cvs:
    # Calculate the results.
    dicc = err_wind_power(str(sys.argv[1]),winds, cv)
    # Fill the dictionaries.
    power['d'].append(dicc['d'])
    power['dplus'].append(dicc['dplus'])

# Save each observed result as a gzipped csv.
results_prefix = './seq_errors/fpr/{0}_'.format(str(sys.argv[1]))

np.savetxt(
    results_prefix+'d.csv.gz',
    [np.array(power['d'])],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)
np.savetxt(
    results_prefix+'dplus.csv.gz',
    [np.array(power['dplus'])],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)
