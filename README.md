# kmer_analysis: evaluation of k-mer based features

## Overview
This project is designed to process proteomic data by splitting the aminoacid-sequences into k-mers and comparing the k-mer frequencies across all input species. It consists of four Python scripts that automate data processing:
- **parse_input.py**: parses the input sequences per species from protein FASTA files
- **kmer_counter.py**: extracts the k-mers from the parsed sequences and counts them
- **analysis.py**: performs exploratory data analysis on the k-mer counts
- **chi2.py**: performs a chi-square test and calculates Cramer's V

The python scripts are embedded into a Snakemake pipeline for easier execution. The pipeline can be run on custom data and the user can choose which k-mer sizes should be generated/included in the analysis.

## Installation
### Prerequisites
- **Python 3.x**
- **Snakemake**
- **Required Python libraries/packages** listed in `requirements.txt`

In case you are working with a Windows operating system instead of Linux, using WSL (Windows Subsystem for Linux) is recommended.

### Setup
1. clone the repository in the desired directory:
   ```bash
   git clone https://github.com/ewearpel/k-mer
   ``` 
2. install dependencies:
   ```bash
   pip install -r ./requirements.txt
   ```
3. install snakemake:
    ```
    pip install snakemake
    ```
    or if using conda:
    ```
    conda install -c bioconda -c conda-forge snakemake
    ```
    Otherwise, check out the installation guide on the official snakemake website: https://snakemake.readthedocs.io/en/v5.6.0/getting_started/installation.html


## Directory Structure
The structure of our directory `k-mer` looks like:
- `bin`: contains the python scripts
- `data`: contains the data to be analyzed as .faa.gz files
- `developer`: contains files for the developer-level documentation
- `README.md`: contains the user-level documentation
- `Snakefile`: contains necessities for the Snakemake pipeline
- `config.yaml`: contains the customizable parts for the Snakemake pipeline
- `requirements.txt`: contains a list of required libraries/packages to execute the code

When cloned from the repository, the folder `data` contains test cases which have to be replaced by the data you want to analyze.

## Usage
### Running the Code
To run the Snakemake pipeline, run the following command from your working directory:
```bash
snakemake -c threads -s /path/to/Snakefile --configfile /path/to/config.yaml 
```
For example, within the working directory the command could look like this:
```bash
snakemake -c 4 --configfile config.yaml
```

If desired, the python-scripts can also be executed separately. Keep in mind to use the correct command line arguments.
```bash
python ./bin/scriptname.py argument1 argument2
```
### Customization
The `config.yaml` file contains the customizable parts of the Snakemake pipeline. Here you can:
- change the **input directory** containing the samples
- change the **samples** to be analyzed (as long as the files are inside the input_dir)
- change the **k-mer size**
- change the **table inclusion threshold** for the chi-square test (recommended to use a number between 1 and 20)

The configuration could look like this:

```
input_dir: data
samples:
  - C.turicensis_non-halo.faa.gz
  - S.aureus_halo.faa.gz
k_values:
  - 3
  - 4
  - 5
table_inclusion_threshold: 5
```

For **analyzing you own data** be sure to:
- add your input file (FASTA and gzipped FASTA format, acquired at e.g. UniProt) to the folder `data` (or an input_dir of your choice)
- change input_dir in the `config.yaml` if necessary
- add the filename to the samples section of `config.yaml`
  
The tool is equipped to handle a multitude of files.

### Results
After running the Snakemake pipeline, the results can be found within the separately created folder `results`. Intermediate results (like the k-mer counts, descriptive statistics, boxplots etc.) can be found within the separately created folder `intermediate`.

**Attention**: When running the Snakemake pipeline, results from previous runs are overwritten, if `config.yaml` has been changed in-between runs. So, if you run it multiple times with different data, be sure to save the results somewhere else. Preferably, you could run the program from a project-specific directory to bypass this.
