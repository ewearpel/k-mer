#!/usr/bin/env python3

import sys
from collections import Counter
import json
import pandas as pd
import os

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

def kmer_counter_multi_input(input_list, k_values):
    """
    Counts k-mers for multiple input files and generates k-mer count contingency tables.

    Args:
        input_list (list of lists): containing one list per input file which in turn contain all sequences of input file
        k_values (list): list of k-mer lengths
        table_inclusion_threshold (int): minimal k-mer frequency per cell for column to be included in contingency table

    Returns:
         result (dict): dictionary containing one contingency table as pd.DataFrame per input k-mer length
    """
    result = {}

    for k in k_values:
        k = int(k) # convert k string from command line to integer

        input_counts_list = [] # set up list to store Counter() objects for each input file
        filenames = [] # store filenames in a list

        for entry in input_list: # iterate over elements in dict representing input files
            counts = Counter() # set up Counter() object
            for seq in entry["sequences"]: # iterate over sequences
                counts.update(kmer_counter(seq, k)) # count k-mers with k_mer counter
            input_counts_list.append(counts) # append Counter() objects to list
            filenames.append(entry["filename"]) # append filenames to list

        # create pd.DataFrame from Counter() objects
        kmer_table = pd.DataFrame.from_records(input_counts_list).fillna(0).astype(int)

        # set index to filenames so that we know which data belongs to which input
        kmer_table.index = filenames


        result[f'k{k}'] = kmer_table # create dictionary of tables

    return result


def count_sequence_length(sequences_by_species):
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

    return sequence_lengths

if __name__ == '__main__':
    seqs_json = sys.argv[1]
    output_dir_k = sys.argv[2]
    output_dir_seq = sys.argv[3]
    k_values = list(sys.argv[4].split(','))

    with open(seqs_json, 'r') as file:
        seqs = json.load(file)

    kmer_counts = kmer_counter_multi_input(seqs, k_values)

    # Save k-mer count tables as TSV files
    for k, df in kmer_counts.items():
        output_file = os.path.join(output_dir_k, f'kmer_counts_{k}.tsv')
        df.to_csv(output_file, sep='\t')
        print(f"{k[1]}-mer count table saved to {output_file}", flush=True)

    print("k-mer count tables have been successfully saved to the folder 'kmer_results'", flush=True)

    # Calculate and save sequence lengths
    sequence_length = count_sequence_length(seqs)

    sequence_length_file = os.path.join(output_dir_seq,'sequence_lengths.json')
    with open(sequence_length_file, 'w') as output:
        json.dump(sequence_length, output)

    print("Sequence lengths have been successfully saved to 'sequence_lengths.json'", flush=True)