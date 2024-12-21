import sys
import json

import numpy as np


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
        if len(species_data) < 2:
            print(f"Error: Expected at least two species for '{kmer_size}', but found {len(species_data)}.")
            sys.exit(1)

        # Parse data for each species
        parsed_data[kmer_size] = {
            "species_1": species_data[0],
            "species_2": species_data[1]
        }

    return parsed_data



def calculate_descriptive_statistics(kmer_data):
    stats = {}
    for kmer_size, species_data in kmer_data.items():
        stats[kmer_size] = {}
        for i, species in enumerate(species_data):
            species_name = f"species_{i + 1}"
            frequencies = list(species.values())
            stats[kmer_size][species_name] = {
                "mean": np.mean(frequencies),
                "median": np.median(frequencies),
                "variance": np.var(frequencies),
                "std_dev": np.std(frequencies),
                "range": (min(frequencies), max(frequencies))
            }
    return stats












"""
if __name__ == "__main__":

    # Load the file path from command-line arguments
    file_path = sys.argv[1]

    # Load and parse the data
    kmer_data = load_and_parse_kmer_data(file_path)

    #save the parsed data
    with open('data_for_analysis.json', 'w') as output:
        json.dump(kmer_data, output)

    # Print the loaded data structure for debugging purposes
    print("Successfully loaded and parsed k-mer data!")
"""