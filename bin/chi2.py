#!/usr/bin/env python3

import sys
from scipy.stats import chi2_contingency
import pandas as pd
import re
import numpy as np
import os

def cramers_v(chi2, n, r, c):
    return np.sqrt(chi2 / (n * min(r - 1, c - 1)))

def chi2_test(tsv_paths, table_inclusion_threshold):
    # store input dataframes in a dictionary
    chi2_dict = {}
    table_inclusion_threshold = int(table_inclusion_threshold)

    for path in tsv_paths:
        # retrieve keys (k-mer length) from file names
        match = re.search(r'k(\d+)', path)
        if match:
            key = f"k{match.group(1)}"

        #set up dictionary to be used to store chi2 results
        chi2_results = {}

        #read tsv files, convert df to dict and store in df_dict with appropriate key
        with open(path, 'r') as file:
            df = pd.read_csv(file, sep='\t')

            numeric_data = df.select_dtypes(include=["number"])
            numeric_data = numeric_data.loc[:, (numeric_data >= int(table_inclusion_threshold)).any(axis=0)]

            chi2, p, dof, expected = chi2_contingency(numeric_data)

            # also calculate cramers v to illustrate association between rows and columns since
            # p-values tend to be 0.0
            n = numeric_data.values.sum()
            r, c = numeric_data.shape
            cramers_v_value = cramers_v(chi2, n, r, c)

            chi2_results['chi2'] = round(float(chi2), 3)
            if p < 1e-2 and p != 0.0:
                chi2_results['p'] = "{:.2e}".format(p)  # wissenschaftliche Notation für sehr kleine p-Werte
            else:
                chi2_results['p'] = round(float(p), 2)
            chi2_results['dof'] = dof
            chi2_results['cramers v'] = round(float(cramers_v_value), 3)

            chi2_dict[key] = chi2_results

    return chi2_dict

def chi2_results_table(chi2_dict, output_path):
    labels = ['k-mer length']
    zero_key = list(chi2_dict.keys())[0]
    add_labels = list(chi2_dict[zero_key].keys())
    labels = labels + add_labels

    with open(output_path, 'w') as file:
        file.write("{:<15} {:<15} {:<15} {:<15} {:<15}\n".format(*labels))

        for key, values in chi2_dict.items():
            row = [key[1]] + list(values.values())
            file.write("{:<15} {:<15} {:<15} {:<15} {:<15}\n".format(*row))


if __name__ == '__main__':
    table_inclusion_threshold = sys.argv[1]
    output_dir = sys.argv[2]
    tsv_paths = sys.argv[3:]


    chi2_dict = chi2_test(tsv_paths,table_inclusion_threshold)
    output_path = os.path.join(output_dir,'chi2_results.txt')
    chi2_results_table(chi2_dict, output_path)
