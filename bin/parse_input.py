#!/usr/bin/env python3

import sys
import gzip
import json
import os

def parse_input(file_paths):
    """
    Reads input files (FASTA format), parses them, and stores contained sequences in a list of lists.

    Args:
        file_paths (list of str): File paths obtained from command line arguments.

    Returns:
        input_data (list of dictionaries): Each dictionary corresponds to input file, contains filename
        and sequences.
    """
    input_data = [] # initializing list to be returned by the function

    for file in file_paths:
        input_seq = [] # initializing inner list representing input file
        current_seq = "" # initializing current sequence
        filename = '.'.join(os.path.splitext(os.path.basename(file))[0].split('.')[:2])

        try:
            # try opening file either as gzip or as plain text file
            with (gzip.open(file, 'rt') if file.endswith(".gz") else open(file, 'r')) as f:
                for line in f: # iterate over lines in file
                    line = line.strip() # cleans lines from spaces, tabs and newline characters
                    if line.startswith(">"): # identify fasta headers
                        if current_seq:
                            input_seq.append(current_seq) # store current sequence to list of sequences
                        current_seq = "" # reset current sequence
                    else:
                        current_seq += line # if no fasta header in the line is detected, append current line to current sequence

                if current_seq:
                    input_seq.append(current_seq) # if no more lines appear, add current sequence to list of sequences

            input_data.append({"filename": filename, "sequences": input_seq})

        # if file cannot be read for any reason
        except (OSError, IOError) as e:
            print(f"Error reading file {file}: {e}")
            input_data.append([])

    return input_data

if __name__ == "__main__":
    # usage info if wrong number of arguments is passed
    if len(sys.argv) < 3: # at least two input files are needed for the pipeline to make any sense
        print("Usage: python parse_input.py <file1> <file2> ... <fileN>")
        sys.exit(1)

    output_file = sys.argv[1]
    filenames = sys.argv[2:]

    print(f"Processing files: {filenames}")
    input_seqs = parse_input(filenames)

    with open(output_file, 'w') as output:
        json.dump(input_seqs, output)  # Save directly to the specified output file

    print(f"Sequences have been successfully saved to '{output_file}'", flush=True)
