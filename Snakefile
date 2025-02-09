configfile: "config.yaml"

rule all:
    input:
        "results/chi2_results.txt"

rule parse_input:
    input:
        expand("{input_dir}/{sample}",sample=config["samples"],input_dir=config["input_dir"])
    output:
        "intermediate/input_seqs.json"
    shell:
        '''
        python3 bin/parse_input.py {output} {input}
        '''

rule count_kmers:
    input:
        rules.parse_input.output
    output:
        tsv_files=expand("intermediate/kmer_results/kmer_counts_k{k}.tsv", k=config["k_values"]),
        sequence_lengths="intermediate/sequence_lengths.json"
    params:
        k_values=lambda wildcards: ','.join(map(str,config["k_values"]))
    shell:
        '''
        python3 bin/kmer_counter.py {input} "intermediate/kmer_results" "intermediate" {params.k_values}
        '''

rule chi2:
    input:
        rules.count_kmers.output.tsv_files
    output:
        "results/chi2_results.txt"
    params:
        threshold = config["table_inclusion_threshold"]
    shell:
        '''
        python3 bin/chi2.py {params.threshold} "results" {input}
        '''