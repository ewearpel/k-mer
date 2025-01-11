import sys
import gzip
import json

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

if __name__ == "__main__":
    filenames = list(sys.argv[1].split(','))

    input_seqs = parse_input(filenames)

    with open('input_seqs.json', 'w') as output:
        json.dump(input_seqs, output)

    print("Sequences have been successfully saved to 'input_seqs.json'")

