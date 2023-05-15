### Dependencies ###
import numpy as np
import sys
### sys.argv[1] = replicate ID ###
### sys.argv[2] = error rate ###

# Define a function to induce sequencing errors.
def sim_seq_errors_windows(error_rate, seq_len=3_000_000_000):
    # Determine the number of sequencing errors.
    seq_errors = round(error_rate * seq_len)
    # Determine the sites where errors will occur.
    p1_all_error_sites = np.unique(np.random.randint(0, seq_len, seq_errors))
    p2_all_error_sites = np.unique(np.random.randint(0, seq_len, seq_errors))
    # Concatenate the error sites for P1 and P2.
    p1_p2_all_error_sites = np.union1d(p1_all_error_sites, p2_all_error_sites)
    # Intialize empty lists to store the results.
    p1_counts = []
    p2_counts = []
    p1_p2_counts = []
    # For every window...
    for wind in range(0, seq_len, 50_000):
        # Determine the number of errors in that window.
        p1_counts.append(np.where(((wind <= p1_all_error_sites) & (p1_all_error_sites < (wind+50_000))))[0].size)
        p2_counts.append(np.where(((wind <= p2_all_error_sites) & (p2_all_error_sites < (wind+50_000))))[0].size)
        p1_p2_counts.append(np.where(((wind <= p1_p2_all_error_sites) & (p1_p2_all_error_sites < (wind+50_000))))[0].size)
    return np.array(p1_counts), np.array(p2_counts), np.array(p1_p2_counts)

# Parse comand-line arguments.
rep_id = int(sys.argv[1])
err_rate = float(sys.argv[2])

# Run the sequencing error simulations.
p1_errors, p2_errors, p1_p2_errors = sim_seq_errors_windows(error_rate=err_rate, seq_len=3_000_000_000)

# Save each observed result as a gzipped csv.
results_prefix = './seq_errors/synthetic_errors/{0}/rep_id_{1}_'.format(err_rate, rep_id)

np.savetxt(
    results_prefix+'p1.csv.gz',
    [p1_errors],
    fmt='%d',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'p2.csv.gz',
    [p2_errors],
    fmt='%d',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'p1_p2.csv.gz',
    [p1_p2_errors],
    fmt='%d',
    delimiter=',',
    newline='\n',
)
