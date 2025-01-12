#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/* ------------------------------------------------------*
 *                                                       *
 *                                                       *
 *              k--mer - Feature Importance              *
 *                                                       *
 *                                                       *
 * ------------------------------------------------------*/


params.input = ""
params.kmer = ""
params.threshold = 1

/*
 * parsing input files
 */
process parseInput {
	input:
	val params.input

	output:
	path "input_seqs.json", emit: input_seqs

	"""
	parse_input.py $params.input
	""" 
}

/*
 * counting k-mers in the input sequences
 */
process kmerCounter {
	input:
	path input_seqs
	val kmer
	val threshold

	output:
	path "kmer_results/*.tsv", emit: kmer_counts
	path "sequence_lengths.json", emit: sequence_lengths

	"""
	echo "DEBUG: kmerCounter received -> $input_seqs"
        echo "DEBUG: k-mer sizes -> $kmer"
        echo "DEBUG: threshold -> $threshold"
	
	kmer_counter.py $input_seqs $params.kmer $params.threshold
	"""
}

workflow {
        if (!params.input) {
                error "No input file provided. Use '--input <input_file>' to specify the input."
        }

        parseInput(params.input)
	
	parseInput.out.input_seqs.view { it -> "DEBUG: parseInput.out.input_seqs -> ${it}" }

	kmerCounter(parseInput.out.input_seqs, params.kmer, params.threshold)
}
