#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/* ------------------------------------------------------*
 *                                                       *
 *                                                       *
 *              k--mer - Feature Importance              *
 *                                                       *
 *                                                       *
 * ------------------------------------------------------*/


params.input = ""       // default to an empty string
params.kmer = ""        // default to an empty string
params.threshold = 1    // default to 1

/*
 * parsing input files
 */
process parseInput {
	input:
	val(input_files)    // uses value of input_files as input for the process

	output:
	path "input_seqs.json", emit: input_seqs    // process produces an output file called "input_seqs.json", labeled as input_seqs for further usage in other processes

	//execute parse_input.py:
	"""
	echo "DEBUG: Received params.input path -> ${params.input}"
	parse_input.py ${input_files.join(' ')}
	"""
}

/*
 * counting k-mers in the input sequences
 */
process kmerCounter {
	input:
	path input_seqs     // for using input_seqs from the previous process
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

process analysis {
    input:
    path kmer_counts
    path sequence_lengths

    output:
    path "descriptive_statistics.json", emit: descriptive_statistics
    path "normalized_descriptive_statistics.json", emit: normalized_descriptive_statistics

    """
    echo "DEBUG: Received TSV paths -> ${kmer_counts.join(' ')}"
        echo "DEBUG: Received sequence lengths path -> $sequence_lengths"

    analysis.py "${kmer_counts.join(' ')}" $sequence_lengths
    """

}

process chi2 {
    input:
    path kmer_counts
    path sequence_lengths

    output:
    path "chi2_results.txt", emit: chi2_results

    """
    chi2.py "${kmer_counts.join(' ')}"
    """
}


workflow {
        if (!params.input) {
                error "No input file provided. Use '--input <input_file>' to specify the input."
        }   // check if an input is provided

	def input_files = params.input.split(' ')	

        parseInput(input_files)

	parseInput.out.input_seqs.view { it -> "DEBUG: parseInput.out.input_seqs -> ${it}" }

	kmerCounter(parseInput.out.input_seqs, params.kmer, params.threshold)

	kmerCounter.out.kmer_counts.view { it -> "DEBUG: kmer_counts -> ${it}" }
	kmerCounter.out.sequence_lengths.view { it -> "DEBUG: kmer_counts -> ${it}" }

	analysis(kmerCounter.out.kmer_counts, kmerCounter.out.sequence_lengths)

	chi2(kmerCounter.out.kmer_counts, kmerCounter.out.sequence_lengths)
}
