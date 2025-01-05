import sys
from collections import Counter
import json
import pandas as pd

def kmer_counter(sequence, k):
    kmer_count = Counter()
    k = int(k)
    for i in range((len(sequence)) - k + 1):
        kmer = sequence[i:i+k]
        kmer_count[kmer] += 1

    return kmer_count

def kmer_counter_multi_input(input_list, k_values, table_inclusion_threshold):
    result = {}
    for k in k_values:
        k = int(k)
        all_inputs_count = []

        input_counts_list = []
        for input in input_list:
            input_counts = Counter()
            for sequence in input:
                input_counts.update(kmer_counter(sequence, k))
            input_counts = dict(sorted(input_counts.items(), key=lambda item: item[1], reverse=True))
            input_counts_list.append(input_counts)

        kmer_table = pd.DataFrame([
            {kmer: counts.get(kmer, 0) for kmer in set().union(*[counts.keys() for counts in input_counts_list])}
            for counts in input_counts_list
        ])
        kmer_table.index = [f"Input_{i + 1}" for i in range(len(input_list))]

        kmer_table = kmer_table.loc[:, (kmer_table >= int(table_inclusion_threshold)).all(axis=0)]

        result[f'k={k}'] = kmer_table
    return result

if __name__ == '__main__':
    seqs_json = sys.argv[1]
    k_values = list(sys.argv[2].split(','))
    table_inclusion_threshold = sys.argv[3]

    with open(seqs_json, 'r') as file:
        seqs = json.load(file)

    kmer_counts = kmer_counter_multi_input(seqs, k_values, table_inclusion_threshold)

    for k, df in kmer_counts.items():
        output_file = f'kmer_counts_{k}.tsv'
        df.to_csv(output_file, sep='\t')
        print(f"{k}-mer count table saved to {output_file}")

    print("k-mer count dictionaries have been successfully saved to 'kmer_counts_k.tsv files.'")