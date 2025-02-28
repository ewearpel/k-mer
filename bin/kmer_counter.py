#!/usr/bin/env python3

"""
kmer_counter.py

This script focuses on the efficient counting of k-mers from biological sequences and provides functionality for generating k-mer count dictionaries.
It supports both single-sequence and multi-sequence processing.

Functions:
    count_kmers(sequence: str, k: int) -> dict:
        Counts k-mers of length `k` in a given biological sequence and returns a dictionary where keys are k-mers and values are their counts.

    merge_kmer_counts(kmer_dicts: list[dict]) -> dict:
        Merges multiple k-mer count dictionaries into a single dictionary, summing the counts for any overlapping k-mers.

    write_kmer_counts_to_file(kmer_counts: dict, output_file: str):
        Writes the k-mer counts to a file in TSV format, with one k-mer and its count per line.

    process_fasta_file(fasta_path: str, k: int) -> dict:
        Processes a FASTA file, counts k-mers for each sequence, and returns a dictionary of k-mer counts.

    kmer_counter_from_fasta_paths(fasta_paths: list[str], k: int, output_dir: str):
        Iterates over multiple FASTA files, computes k-mer counts for each, and saves them in TSV format.

Main Execution:
    Command-line arguments:
        1. fasta_paths (list of str): List of paths to the FASTA files to be processed.
        2. k (int): Length of k-mers to be counted.
        3. output_dir (str): Directory where k-mer count TSV files will be saved.

    Workflow:
        1. Read biological sequences from FASTA files.
        2. Count k-mers for each sequence.
        3. Save the results to the specified output directory in TSV format.
"""

import sys
import json
import pandas as pd
import os
import numpy as np
from collections import Counter
from multiprocessing import Pool, cpu_count

def kmer_counter(sequence, k):
    """
    Counts k-mers of defined length and stores counts in a dictionary.

    Args:
        input sequence (str)
        k (str): k-mer length

    Returns:
        kmer_count (dict): {k-mer: count}
    """
    kmer_count = Counter() # sets up kmer_count as Counter object which is a subclass of dict
    k = int(k) # converts input string for k-mer length to integer
    for i in range((len(sequence)) - k + 1): # iterates over the entire sequence
        kmer = sequence[i:i+k] # defines kmer per iteration
        kmer_count[kmer] += 1 # value belonging to kmer key is incremented

    return kmer_count


def kmer_counter_multi_input(input_list, k_values, output_dir):
    """
    Counts k-mers for multiple input files and generates k-mer count contingency tables.

    Args:
        input_list (list of lists): containing one list per input file which in turn contain all sequences of input file
        k_values (list): list of k-mer lengths
        table_inclusion_threshold (int): minimal k-mer frequency per cell for column to be included in contingency table

    Returns:
         .tsv files containing contingency table for each k-mer length
    """
    for k in k_values:
        k = int(k) # convert k string from command line to integer

        input_counts_list = [] # set up list to store Counter() objects for each input file
        filenames = [] # store filenames in a list

        for entry in input_list: # iterate over elements in dict representing input files
            counts = Counter(kmer for seq in entry["sequences"] for kmer, count in kmer_counter(seq, k).items() for _ in range(count))
            input_counts_list.append(counts) # append Counter() objects to list
            filenames.append(entry["filename"]) # append filenames to list

        # create pd.DataFrame from Counter() objects and save to tsv file
        kmer_table = pd.DataFrame.from_records(input_counts_list).fillna(0).astype(int)
        kmer_table.index = filenames
        output_file = os.path.join(output_dir, f'kmer_counts_k{k}.tsv')
        kmer_table.to_csv(output_file, sep='\t')
        print(f"{k}-mer count table saved to {output_file}", flush=True)



def count_sequence_length(sequences_by_species, output_dir_seq):
    """
        Calculates the total sequence length for each input file.

        Args:
            sequences_by_species (list): List of dictionaries with "filename" and "sequences" keys.

        Returns:
            sequence_lengths (dict): Dictionary with filenames as keys and total sequence lengths as values.
        """
    sequence_lengths = {
        entry["filename"]: sum(len(sequence) for sequence in entry["sequences"])
        for entry in sequences_by_species
    }

    sequence_length_file = os.path.join(output_dir_seq, 'sequence_lengths.json')
    with open(sequence_length_file, 'w') as output:
        json.dump(sequence_lengths, output)

    print("Sequence lengths saved to 'sequence_lengths.json'", flush=True)

if __name__ == '__main__':
    seqs_json = sys.argv[1]
    output_dir_k = sys.argv[2]
    output_dir_seq = sys.argv[3]
    k_values = list(sys.argv[4].split(','))

    if len(sys.argv) < 5:
        print(
            "Usage: python kmer_counter.py <input_sequences.json> <output_dir_k> <output_dir_seq> <k-mer_values_comma_separated>")
        sys.exit(1)

    with open(seqs_json, 'r') as file:
        seqs = json.load(file)

    kmer_counter_multi_input(seqs, k_values, output_dir_k)
    count_sequence_length(seqs, output_dir_seq)