import sys
import json
import numpy as np
import pandas as pd
import os
import glob
from scipy.stats import chi2_contingency

def calculate_descriptive_statistics(kmer_data):
    stats = {}
    normalized_stats = {}

    for tsv_file in tsv_paths:
        kmer_size = os.path.basename(tsv_file).split('_')[2]  # extract k-mer size from file name
        calculation_df = pd.read_csv(tsv_file, sep='\t', index_col=0)

        stats[kmer_size] = {}
        normalized_stats[kmer_size] = {}

        for species_name in calculation_df.index:
            frequencies = calculation_df.loc[species_name].values # Extract the frequencies

            stats[kmer_size][species_name] = {
                "mean": float(np.mean(frequencies)),    # float() because numpy types are not serializable in json
                "median": float(np.median(frequencies)),
                "variance": float(np.var(frequencies)),
                "std_dev": float(np.std(frequencies)),
                "max_count": float(max(frequencies))
            }

            normalized_frequencies = []
            total_length = sequence_lengths.get(species_name, 1) # 1 as default to avoid dividing by zero
            for freq in frequencies:
                normalized_frequencies.append(freq / total_length)

            normalized_stats[kmer_size][species_name] = {
                "mean": float(np.mean(normalized_frequencies)),
                "median": float(np.median(normalized_frequencies)),
                "variance": float(np.var(normalized_frequencies)),
                "std_dev": float(np.std(normalized_frequencies)),
                "max_count": float(max(normalized_frequencies))
            }

    return stats, normalized_stats


def statistics_to_table(calculated_statistics):
    # initialize table structure
    table_data = {
        "Statistic": ["mean", "median", "variance", "std_dev", "max_count"]  # Row names
    }

    # fill the columns dynamically from JSON
    for kmer_size, species_data in calculated_statistics.items():
        for species_name, species_stats in species_data.items():  # Avoid overwriting `stats`
            column_name = f"{kmer_size} for {species_name}"
            table_data[column_name] = [
                species_stats["mean"],
                species_stats["median"],
                species_stats["variance"],
                species_stats["std_dev"],
                species_stats["max_count"]
            ]

    # Create a DataFrame
    table_df = pd.DataFrame(table_data)

    return table_df

"""
# because json.dump cannot serialize numpy type numbers
def convert_numpy_to_python_types(data):
    if isinstance(data, dict):
        return {key: convert_numpy_to_python_types(value) for key, value in data.items()}
    elif isinstance(data, list):
        return [convert_numpy_to_python_types(item) for item in data]
    elif isinstance(data, (np.integer, np.floating)):
        return data.item()  # convert to native Python type
    else:
        return data
"""


"""
def chi_square_test(kmer_data):
    chi_square_results = {}

    for kmer_size, species_data in kmer_data.items():
        # prepare contingency table
        contingency_table = []
        species_names = []
        filtered_kmers = set()

        # collect all k-mers that have non-zero counts across all species
        for species_name, kmer_frequencies in species_data.items():
            for kmer, count in kmer_frequencies.items():
                if count > 0: # add kmer when it's count is above 0
                    filtered_kmers.add(kmer)
        filtered_kmers = sorted(filtered_kmers)

        # build the contingency table for valid k-mers with normalized frequencies
        for species_name, kmer_frequencies in species_data.items():
            normalized_frequencies = []
            total_length = sequence_lengths.get(species_name, 1)

            for kmer in filtered_kmers:
                freq = kmer_frequencies.get(kmer, 0)
                normalized_frequencies.append(freq / total_length)

            contingency_table.append(normalized_frequencies)
            species_names.append(species_name)

        # perform chi-square test
        try:
            chi2, p, dof, expected = chi2_contingency(contingency_table)
            chi_square_results[kmer_size] = {
                "chi2_statistic": chi2,
                "p_value": p,
                "degrees_of_freedom": dof,
                "species_names": species_names
            }
        except ValueError as e:
            chi_square_results[kmer_size] = {
                "error": str(e),
                "species_names": species_names
            }

    return chi_square_results

def chi2_to_table(chi_square_results):
    # initialize table structure
    table_data = {
        "results": ["chi2_statistic", "p_value", "degrees_of_freedom", "species_names"]  # Row names
    }

    for kmer_size, result in chi_square_results.items():
        column_name = f"{kmer_size}"
        table_data[column_name] = [
            result["chi2_statistic"],
            result["p_value"],
            result["degrees_of_freedom"],
            result["species_names"]
        ]

    # Create a DataFrame
    chi2_table_df = pd.DataFrame(table_data)

    return chi2_table_df
"""



if __name__ == "__main__":

    # Load the file path from command-line arguments
    folder_path = sys.argv[1]
    tsv_paths = glob.glob(os.path.join(folder_path, '*.tsv'))
    lengths_path = sys.argv[2]

    if not tsv_paths:
        print(f"No tsv-files found in the folder '{folder_path}'.")
        sys.exit(1)

    with open(lengths_path, 'r') as file:
        sequence_lengths = json.load(file)

    # calculate statistics
    stats, normalized_stats = calculate_descriptive_statistics(tsv_paths)

    with open('descriptive_statistics.json', 'w') as output:
        json.dump(stats, output)

    with open('normalized_descriptive_statistics.json', 'w') as output:
        json.dump(normalized_stats, output)

    print("Successfully calculated descriptive statistics!")

    # put statistics in a table
    absolute_table = statistics_to_table(stats)
    normalized_table = statistics_to_table(normalized_stats)

    # save as csv?
    # absolute_table.to_csv("absolute_statistics_table.csv", index=False)
    # normalized_table.to_csv("normalized_statistics_table.csv", index=False)

    print("\nabsolute statistics:") # \n for better visibility
    print(absolute_table.to_string(index=False))

    print("\nnormalized statistics:")
    print(normalized_table.to_string(index=False))


