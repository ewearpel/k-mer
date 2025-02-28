#!/usr/bin/env python3

# Author: Sophie Wieser

"""
analysis.py

This script calculates descriptive statistics for k-mer counts across species and generates boxplots for visualization.
It processes k-mer count TSV files and sequence lengths to produce both raw and normalized statistics.

Functions:
    calculate_descriptive_statistics(kmer_data) -> tuple[dict, dict]:
        Calculate descriptive statistics (mean, median, variance, std dev, max count) for each species
        from k-mer count data, both raw and normalized by sequence length.

    statistics_to_table(calculated_statistics: dict) -> pandas.DataFrame:
        Convert the calculated statistics into a pandas DataFrame for easy export or reporting.

    plot_statistic_boxplots(calculated_statistics: dict, output_folder: str, title_prefix: str = "Descriptive Statistics"):
        Generate boxplots for each k-mer size across species and save them as PNG files.

Main Execution:
    Command-line arguments:
        1. output_dir_box (str): Directory for saving boxplots.
        2. output_dir_stats (str): Directory for saving statistics in JSON format.
        3. lengths_path (str): JSON file path containing the sequence lengths for each species.
        4. tsv_paths (list of str): Paths to the k-mer count TSV files.

    Workflow:
        1. Load sequence lengths and k-mer counts.
        2. Calculate descriptive statistics.
        3. Save statistics as JSON files.
        4. Plot boxplots for normalized statistics.
"""

import sys
import json
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns

import glob


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

def plot_statistic_boxplots(calculated_statistics, output_folder,  title_prefix="Descriptive Statistics"):


    # Convert the nested dictionary into a long-form DataFrame
    rows = []
    for kmer_size, species_data in calculated_statistics.items():
        for species_name, stats in species_data.items():
            for stat_name, stat_value in stats.items():
                rows.append({
                    "kmer_size": kmer_size,
                    "species_name": species_name,
                    "statistic": stat_name,
                    "value": stat_value
                })

    stats_df = pd.DataFrame(rows)

    # Get all distinct k-mer sizes from the DataFrame
    unique_kmers = stats_df["kmer_size"].unique()

    # For each k-mer size, plot a boxplot of "value" across the different species
    for kmer in unique_kmers:
        subset_df = stats_df[stats_df["kmer_size"] == kmer]

        # If kmer is like "k=2.tsv", remove the ".tsv" part:
        if kmer.endswith(".tsv"):
            kmer_base = kmer.replace(".tsv", "")
        else:
            kmer_base = kmer

        plt.figure(figsize=(6, 12))
        # x="species_name" -> side-by-side boxes for each species
        sns.boxplot(x="species_name", y="value", data=subset_df)
        plt.title(f"{title_prefix} for k-mer size: {kmer_base}")
        plt.xlabel("Species")
        plt.ylabel("Value")

        # Rotate x-labels if you have many species names
        plt.xticks(rotation=45)
        plt.tight_layout()
        plot_filename = os.path.join(output_folder, f"boxplots_{kmer_base}.png")
        plt.savefig(plot_filename)



if __name__ == "__main__":
    # Load the file path from command-line arguments
    output_dir_box = sys.argv[1]
    output_dir_stats = sys.argv[2]
    lengths_path = sys.argv[3]
    tsv_paths = sys.argv[4:]

    if len(sys.argv) < 5:
        print("Usage: python analysis.py <output_dir_box> <output_dir_stats> <lengths_path.json> <kmer_tsv_paths...>")
        sys.exit(1)

    with open(lengths_path, 'r') as file:
        sequence_lengths = json.load(file)

    # calculate statistics
    stats, normalized_stats = calculate_descriptive_statistics(tsv_paths)

    descriptive_path = os.path.join(output_dir_stats, 'descriptive_statistics.json')
    with open(descriptive_path, 'w') as output:
        json.dump(stats, output)

    norm_descriptive_path = os.path.join(output_dir_stats,'normalized_descriptive_statistics.json')
    with open(norm_descriptive_path, 'w') as output:
        json.dump(normalized_stats, output)

    print("Successfully calculated descriptive statistics!")

    # put statistics in a table
    absolute_table = statistics_to_table(stats)
    normalized_table = statistics_to_table(normalized_stats)

    plot_statistic_boxplots(normalized_stats, output_dir_box)




