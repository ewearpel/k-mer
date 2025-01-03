import sys
import json
from scipy.stats import chi2_contingency
import pandas as pd

def filter_kmer_for_chi2()


if __name__ == '__main__':
    kmer_dict_path = sys.argv[1]

    with open(kmer_dict_path, 'r') as file:
        content = file.read()
        kmer_dict = json.loads(content)