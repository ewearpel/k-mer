import sys
import json
from scipy.stats import chi2_contingency
import pandas as pd
import re



def chi2_test(tsv_paths):
    # store input dataframes in a dictionary
    chi2_dict = {}

    for path in tsv_paths:
        # retrieve keys (k-mer length) from file names
        match = re.search(r'k=(\d+)', path)
        if match:
            key = f"k={match.group(1)}"

        #set up dictionary to be used to store chi2 results
        chi2_results = {}

        #read tsv files, convert df to dict and store in df_dict with appropriate key
        with open(path, 'r') as file:
            df = pd.read_csv(file, sep='\t')
            numeric_data = df.select_dtypes(include=["number"])
            print(numeric_data)

            chi2, p, dof, expected = chi2_contingency(numeric_data)
            chi2_results['chi2'] = float(chi2)
            chi2_results['p'] = float(p)
            chi2_results['dof'] = dof
            chi2_results['expected'] = expected

            chi2_dict[key] = chi2_results

    return chi2_dict

if __name__ == '__main__':
    tsv_paths = list(sys.argv[1].split(','))

    chi2_dict = chi2_test(tsv_paths)
    print(chi2_dict)