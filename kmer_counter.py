import copy
import sys
import json
from itertools import product

def all_kmer_combinations(size):
    possible_aminoacids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'Y']
    # include U?
    combinations = [''.join(letter) for letter in product(possible_aminoacids, repeat = size)]
    empty_dictionary = {string: 0 for string in combinations}

    return empty_dictionary

def update_count(counts,sequence,k):
    for i in range((len(sequence)) - k + 1):
        kmer = sequence[i:i + k]
        counts[kmer] += 1

def kmer_counter_multi_input(sequences_by_species,k_values):
    result = {}
    for k in k_values:
        k = int(k)
        all_inputs_count = []
        possible_combinations = all_kmer_combinations(k)
        for species_sequences in sequences_by_species:
            counts_per_species = copy.deepcopy(possible_combinations)
            for sequence in species_sequences:
                update_count(counts_per_species,sequence,k)
            all_inputs_count.append(counts_per_species)

        result[f'k={k}'] = all_inputs_count

    return result

def count_sequence_length(sequences_by_species):
    sequence_lengths = {}
    for index,species in enumerate(sequences_by_species):
        length = 0
        for sequence in species:
            length = length + len(sequence)
        sequence_lengths[f'species_{index+1}'] = length

    print(sequence_lengths)
    return sequence_lengths


if __name__ == '__main__':
    seqs_json = sys.argv[1]
    k_values = list(sys.argv[2].split(','))

    with open(seqs_json, 'r') as file:
        seqs = json.load(file)

    kmer_counts = kmer_counter_multi_input(seqs, k_values)

    with open('kmer_counts.json', 'w') as output:
        json.dump(kmer_counts, output)

    print("k-mer count dictionaries have been successfully saved to 'kmer_counts.json'")

    sequence_length = count_sequence_length(seqs)

    with open('sequence_lengths.json', 'w') as output:
        json.dump(sequence_length, output)

