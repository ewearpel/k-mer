import sys
from collections import Counter
import json

def kmer_counter(sequence, k):
    kmer_count = Counter()
    k = int(k)
    for i in range((len(sequence)) - k + 1):
        kmer = sequence[i:i+k]
        kmer_count[kmer] += 1

    return kmer_count

def kmer_counter_multi_input(input_list,k_values):
    result = {}
    for k in k_values:
        k = int(k)
        all_inputs_count = []
        for input in input_list:
            input_counts = Counter()
            for sequence in input:
                input_counts.update(kmer_counter(sequence, k))
            all_inputs_count.append(dict(input_counts))
        result[f'k={k}'] = all_inputs_count

    return result

if __name__ == '__main__':
    seqs_json = sys.argv[1]
    k_values = list(sys.argv[2].split(','))

    with open(seqs_json, 'r') as file:
        seqs = json.load(file)

    kmer_counts = kmer_counter_multi_input(seqs, k_values)

    with open('kmer_counts.json', 'w') as output:
        json.dump(kmer_counts, output)

    print("k-mer count dictionaries have been successfully saved to 'kmer_counts.json'")