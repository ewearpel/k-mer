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
            print(total_length)
            for freq in frequencies:
                normalized_frequencies.append(freq / total_length)

            normalized_stats[kmer_size][species_name] = {
                "mean": np.mean(normalized_frequencies),
                "median": np.median(normalized_frequencies),
                "variance": np.var(normalized_frequencies),
                "std_dev": np.std(normalized_frequencies),
                "range": (min(normalized_frequencies), max(normalized_frequencies))
            }

    return stats, normalized_stats


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










"""
def calculate_descriptive_statistics(kmer_data):
    stats = {}
    for kmer_size, species_data in kmer_data.items():
        stats[kmer_size] = {}
        for i, species in enumerate(species_data):
            species_name = f"species_{i + 1}"
            print(species_data)
            # frequencies = list(species_data.values())
            frequencies = list(kmer_frequencies.values())
            stats[kmer_size][species_name] = {
                "mean": np.mean(frequencies),
                "median": np.median(frequencies),
                "variance": np.var(frequencies),
                "std_dev": np.std(frequencies),
                "range": (min(frequencies), max(frequencies))
            }
    return stats
"""