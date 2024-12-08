import sys
import gzip
from collections import Counter
import pandas as pd


def parse_input(file_paths):
    input_seqs = []

    for file in file_paths:
        input_seq = []
        current_seq = ''
        file_handle = gzip.open(file, 'rt') if file.endswith(".gz") else open(file, 'r')

        with file_handle as f:
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if current_seq:
                        input_seq.append(current_seq)
                        current_seq = ''
                    else:
                        current_seq = ''
                else:
                    current_seq += line
        input_seqs.append(input_seq)

    return input_seqs

def kmer_counter(sequence, k):
    kmer_count = Counter()
    k = int(k)
    for i in range((len(sequence)) - k + 1):
        kmer = sequence[i:i+k]
        kmer_count[kmer] += 1

    return kmer_count

def kmer_counter_multi_input(input_list,k):
    all_inputs_count = []
    k = int(k)
    for input in input_list:
        input_counts = Counter()
        for sequence in input:
            input_counts.update(kmer_counter(sequence,k))
        all_inputs_count.append(input_counts)

    return all_inputs_count

def create_dataframe(all_inputs_count):
    kmers = set(kmer for counts in all_inputs_count for kmer in counts)
    pass

if __name__ == "__main__":
    filenames = list(sys.argv[1].split(','))
    kmer_lengths = list(sys.argv[2].split(','))

    inputs = parse_input(filenames)
    for k in kmer_lengths:
        counts = kmer_counter_multi_input(inputs,k)

