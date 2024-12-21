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

def kmer_counter_multi_input(sequences_by_species,k_values):
    result = {}
    for k in k_values:
        k = int(k)
        all_inputs_count = []
        for species_sequences in sequences_by_species:
            counts_per_species = Counter()
            for sequence in species_sequences:
                counts_per_species.update(kmer_counter(sequence, k))

            # sort the k-mers by highest count
            dictionary = dict(counts_per_species)
            sorted_dictionary = dict(sorted(dictionary.items(),reverse=True,key=lambda item: item[1]))
            all_inputs_count.append(sorted_dictionary)

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