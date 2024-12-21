# k-mer project
## parse input
To parse FASTA-files containing multiple sequences, you need to use "parse_input.py" and sepcify the input data. ".gz" folders are unzipped automatically. Multiple input files can be added, but need to be separated via a comma.

Example command line input for two FASTA-files (within the correct folder):
Python parse_input.py ./data/C.turicensis_non-halo.faa.gz,./data/S.aureus_halo.faa.gz

The python script extracts the particular sequences and outputs them within the current folder as new file. If the parsing was successful, the statement "Sequences have been successfully saved to 'input_seqs.json'" is printed. This file contains the parsed sequences as:
\[[sequences from C.turicensis], [sequences from S.aureus]]

## extract and count k-mers
To extract and count the k-mers from the newly created sequence-file, you need to use "kmer_counter.py". You need to specify the k-mer length. Multiple lengths need to be separated via a comma.

Example command line input for k-mers of length 2 and 3:
Python kmer_counter.py ./input_seqs.json 2,3

This python script extracts k-mers from the input sequences, counts them per species and outputs the data within the current folder as a new file. If the extraction and counting was successful, the statement "k-mer count dictionaries have been successfully saved to 'kmer_counts.json'" is printed. This file contains the results as:
{"k=2": [{2-mers of C.turicensis},{2-mers of S.aureus}], "k=3": [{3-mers of C.turicensis},{3-mers of S.aureus}]}
{"k=2": species1: {2-mers of C.turicensis}, species2:{2-mers of S.aureus}], "k=3": [{3-mers of C.turicensis},{3-mers of S.aureus}]}

