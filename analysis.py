import sys
import json
import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency


def load_and_parse_kmer_data(file_path):
    try:
        with open(file_path, 'r') as file:
            data = json.load(file)
    # optional except-statements for better user experience
    except FileNotFoundError:
        print(f"Error: File '{file_path}' not found.")
        sys.exit(1) #terminates the program immediately
    except json.JSONDecodeError:
        print(f"Error: File '{file_path}' is not a valid JSON file.")
        sys.exit(1)

    parsed_data = {}
    for kmer_size, species_data in data.items():
        parsed_data[kmer_size] = {
            "species_1": species_data[0],
            "species_2": species_data[1]
        }

    return parsed_data


def calculate_descriptive_statistics(kmer_data):
    stats = {}
    normalized_stats = {}
    for kmer_size, species_data in kmer_data.items():
        stats[kmer_size] = {}
        normalized_stats[kmer_size] = {}
        for species_name, kmer_frequencies in species_data.items():
            # Ensure that kmer_frequencies is a dictionary
            if not isinstance(kmer_frequencies, dict):
                print(f"Error: Unexpected data type for {species_name} in {kmer_size}.")
                continue

            frequencies = list(kmer_frequencies.values()) # Extract the frequencies


            stats[kmer_size][species_name] = {
                "mean": np.mean(frequencies),
                "median": np.median(frequencies),
                "variance": np.var(frequencies),
                "std_dev": np.std(frequencies),
                "range": (min(frequencies), max(frequencies))
            }

            normalized_frequencies = []
            total_length = sequence_lengths.get(species_name, 1) # 1 as default to avoid dividing by zero
            for freq in frequencies:
                normalized_frequencies.append(freq / total_length)

            normalized_stats[kmer_size][species_name] = {
                "mean": np.mean(normalized_frequencies),
                "median": np.median(normalized_frequencies),
                "variance": np.var(normalized_frequencies),
                "std_dev": np.std(normalized_frequencies)
                #, "range": (min(normalized_frequencies), max(normalized_frequencies))
            }

    return stats, normalized_stats


def statistics_to_table(calculated_statistics):
    # initialize table structure
    table_data = {
        "Statistic": ["mean", "median", "variance", "std_dev"]  # Row names
    }

    # fill the columns dynamically from JSON
    for kmer_size, species_data in calculated_statistics.items():
        for species_name, species_stats in species_data.items():  # Avoid overwriting `stats`
            column_name = f"{kmer_size} for {species_name}"
            table_data[column_name] = [
                species_stats["mean"],
                species_stats["median"],
                species_stats["variance"],
                species_stats["std_dev"]
                #, species_stats["range"]
            ]

    # Create a DataFrame
    table_df = pd.DataFrame(table_data)

    return table_df

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



if __name__ == "__main__":

    # Load the file path from command-line arguments
    file_path = sys.argv[1]
    lengths_path = sys.argv[2]

    with open(lengths_path, 'r') as file:
        sequence_lengths = json.load(file)

    # Load and parse the data
    kmer_data = load_and_parse_kmer_data(file_path)

    #save the parsed data
    with open('data_for_analysis.json', 'w') as output:
        json.dump(kmer_data, output)

    print("Successfully loaded and parsed k-mer data!")

    # calculate statistics
    stats, normalized_stats = calculate_descriptive_statistics(kmer_data)

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

    chi_square = chi_square_test(kmer_data)

    chi2_table = chi2_to_table(chi_square)
    print("\nchi square results:")
    print(chi2_table.to_string(index=False))

