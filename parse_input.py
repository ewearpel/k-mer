import sys
import gzip

def parse_input(file_paths):
    proteomes = []

    for file in file_paths:
        proteome = []
        current_seq = ''
        file_handle = gzip.open(file, 'rt') if file.endswith(".gz") else open(file, 'r')

        with file_handle as f:
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if current_seq:
                        proteome.append(current_seq)
                        current_seq = ''
                    else:
                        current_seq = ''
                else:
                    current_seq += line
        proteomes.append(proteome)

    return proteomes

if __name__ == "__main__":
    filenames = list(sys.argv[1].split(','))

    proteomes=parse_input(filenames)
    print(proteomes)