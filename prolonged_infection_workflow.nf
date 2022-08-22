#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// NOTE that helpful runtime parameters are in the file `nextflow.config`
// We recommend that you modify any file paths, inputs, or outputs there.


// Defining workflow for generating non-supplementary figures for Halfmann et al. 2022
workflow {

	// INPUT CHANNELS
	ch_consensus_seqs = Channel
		.fromPath( params.consensus_seqs )

	ch_ont_reads = Channel // this channel streams in data on the timepoints sequenced on ONT instruments
		.fromPath ( params.primer_key )
		.splitCsv ( sep: ",", header: true )
		.map      { row -> tuple(row.sra_id, row.file_basename, row.platform, file(row.primer_set)) }
		.filter   { it[2] == "ont" }

	ch_illumina_reads = Channel // this channel streams in data on the timepoints sequenced on Illumina instruments
		.fromPath ( params.primer_key )
		.splitCsv ( sep: ",", header: true )
		.map      { row -> tuple(row.sra_id, row.file_basename, row.platform, file(row.primer_set)) }
		.filter   { it[2] == "illumina" }


	// WORKFLOW STEPS
	// Classifying pango lineages
	PANGO_CLASSIFICATION (
		ch_consensus_seqs
	)

	// Plotting Ct values, figure 1A
	FIGURE_1A_PLOTTING ()

	// Processing & variant-calling consensus sequences for figure 2B
	ISOLATING_CONSENSUS_SEQS (
		ch_consensus_seqs
	)

	CONSENSUS_ALIGNMENT (
		ISOLATING_CONSENSUS_SEQS.out
			.flatten()
			.map { file -> tuple(file.simpleName, file) }
	)

	CONSENSUS_PILEUP (
		CONSENSUS_ALIGNMENT.out
	)

	CONSENSUS_VARIANT_CALLING (
		CONSENSUS_PILEUP.out
	)

	// Processing reads for defining amplicon dropouts, which will also be plotted in figure 2B
	// ONT read portion
	GET_ONT_READS (
		ch_ont_reads
	)

	ONT_READ_MAPPING (
		GET_ONT_READS.out
	)

	ONT_ALIGNMENT_POLISHING (
		ONT_READ_MAPPING.out
	)

	ONT_LOWCOV_ANNOTATION (
		ONT_ALIGNMENT_POLISHING.out.bam
	)

	// Illumina read portion
	GET_ILL_READS (
		ch_illumina_reads
	)

	ILL_READ_MAPPING (
		GET_ILL_READS.out
	)

	ILL_ALIGNMENT_POLISHING (
		ILL_READ_MAPPING.out
	)

	ILL_LOWCOV_ANNOTATION (
		ILL_ALIGNMENT_POLISHING.out.bam
	)

	// plotting figure 2A
	FIGURE_2A_PLOTTING (
		ch_consensus_seqs,
		ONT_LOWCOV_ANNOTATION.out.collect(),
		ILL_LOWCOV_ANNOTATION.out, // we would use collect() here, but there's only one Illumina timepoint
		CONSENSUS_VARIANT_CALLING.out.collect()
	)
	
	TAR_FIG2A_DATA (
		FIGURE_2A_PLOTTING.out.cue
	)

	// Plotting neutralization assay results, figure 2C
	FIGURE_2C_PLOTTING ()

}


process PANGO_CLASSIFICATION {

	// Classifying pango lineages for each timepoint

	publishDir params.results, mode: 'move'

	input:
	file(consensus)

	output:
	file("prolonged_infection_lineage_report.csv")

	script:
	"""
	pangolin --outfile 'prolonged_infection_lineage_report.csv' ${consensus}
	"""
}


process FIGURE_1A_PLOTTING {

	// Plotting Ct values

	publishDir params.visuals, pattern: '*.pdf', mode: 'move'

	output:
	file("*.pdf")

	script:
	"""
	Figure1A_Ct_through_time.R ${params.Ct_data}
	"""

}


process ISOLATING_CONSENSUS_SEQS {

	// Taking each consensus sequence into its own, separate FASTA

	input:
	path(consensus)

	output:
	file("USA*.fasta")

	script:
	"""
	fasta_sep.R ${consensus}
	"""

}


process CONSENSUS_ALIGNMENT {

	// Aligning consensus sequences so that they are in sam format

	tag "${timepoint}"

	input:
	tuple val(timepoint), file(fasta)

	output:
	tuple val(timepoint), path("*.sam")

	script:
	"""
	minimap2 -a ${params.refseq} ${fasta} > ${timepoint}.sam
	"""

}


process CONSENSUS_PILEUP {

	// Creating pileups to be piped into iVar in the next process

	tag "${timepoint}"

	input:
	tuple val(timepoint), path(sam)

	output:
	tuple val(timepoint), path("*.mpileup")

	script:
	"""
	cat ${sam} \
	  | samtools view -Sb - \
	  | samtools sort - > tempfile
	  samtools mpileup -aa -f ${params.refseq} --output ${timepoint}.mpileup tempfile
	"""

}


process CONSENSUS_VARIANT_CALLING {

	// Calling variants and protein effects for those variants with iVar

	tag "${timepoint}"
	
	publishDir params.results_data_files, pattern: '*_consensus_variant_table.tsv', mode: 'copy'

	input:
	tuple val(timepoint), path(mpileup)

	output:
	path("*_consensus_variant_table.tsv")

	script:
	"""
	cat ${mpileup} \
	  | ivar variants -p ${timepoint}_consensus_variant_table \
	  -t 0 -m 1 -q 1 -r ${params.refseq} -g ${params.refgff}
	"""

}


process GET_ONT_READS {

	// Using sra-tools fasterq-dump to full these files now, though note:
	// https://github.com/ncbi/sra-tools/issues/463

	tag "${timepoint}"
	
	publishDir params.results_data_files, pattern: '*.fastq.gz', mode: 'copy'
	
	errorStrategy 'retry'
	maxRetries 2

	input:
	tuple val(sra_id), val(timepoint), val(platform), file(primers)

	output:
	tuple val(timepoint), file("*.fastq.gz"), file(primers)

	script:
	"""

	prefetch ${sra_id}
	fasterq-dump ${sra_id}/${sra_id}.sra \
	--concatenate-reads --skip-technical --quiet && \
	gzip ${sra_id}.sra.fastq
	mv ${sra_id}.sra.fastq.gz ${timepoint}.fastq.gz
	rm -rf ${sra_id}/
	rm -rf fasterq.tmp.*
	
	"""

}


process ONT_READ_MAPPING {

	// Mapping ONT reads from each timepoint

	tag "${timepoint}"
	
	memory { 2.GB * task.attempt }
	errorStrategy 'retry'
	maxRetries 2
	// containerOptions '--memory=2g'
	cpus 1

	input:
	tuple val(timepoint), file(fastq), file(primers)

	output:
	tuple val(timepoint), file("*.sam"), file(primers)

	script:
	"""
	minimap2 -ax map-ont -t 1 -o ${timepoint}.sam ${params.refseq} ${fastq} && \
	rm -f `realpath *.fastq.gz`
	"""

}


process ONT_ALIGNMENT_POLISHING {

	// Using SAMTOOLS to clean and prep alignments

	tag "${timepoint}"
	
	publishDir params.results_data_files, pattern: '*.bam', mode: 'copy'
	publishDir params.results_data_files, pattern: '*.bam.bai', mode: 'move'

	input:
	tuple val(timepoint), file(sam), file(primers)

	output:
	tuple val(timepoint), file("*.bam"), emit: bam
	path("*.bam.bai"), emit: bai

	script:
	"""

	cat ${sam} \
	| samtools ampliconclip -b ${primers} - \
	| samtools view -Sb - \
	| samtools sort - > ${timepoint}.bam && \
	rm -f `realpath ${sam}`

	samtools index ${timepoint}.bam

	"""

}

process ONT_LOWCOV_ANNOTATION {

	// Creating BED files that annotate regions of each genome with less than 20x coverage

	tag "${timepoint}"
	
	publishDir params.results_data_files, mode: 'copy'
	
	input:
	tuple val(timepoint), file(bam)

	output:
	file("*.bed")

	script:
	"""
	covtobed --max-cov=20 ${bam} > ${timepoint}.bed && \
	rm -f `realpath *.bam`
	"""

}


process GET_ILL_READS {

	tag "${timepoint}"
	
	publishDir params.results_data_files, pattern: '*.fastq.gz', mode: 'copy'

	input:
	tuple val(sra_id), val(timepoint), val(platform), file(primers)

	output:
	tuple val(timepoint), file("*_R1.fastq.gz"), file("*_R2.fastq.gz"), file(primers)

	script:
	"""

	prefetch ${sra_id}
	fasterq-dump ${sra_id}/${sra_id}.sra \
	--split-files --skip-technical --quiet && \
	gzip ${sra_id}.sra_1.fastq && gzip ${sra_id}.sra_2.fastq
	mv ${sra_id}.sra_1.fastq.gz ${timepoint}_R1.fastq.gz
	mv ${sra_id}.sra_2.fastq.gz ${timepoint}_R2.fastq.gz
	rm -rf ${sra_id}/
	rm -rf fasterq.tmp.*

	"""

}


process ILL_READ_MAPPING {

	// Mapping ONT reads from each timepoint

	tag "${timepoint}"
	
	memory { 2.GB * task.attempt }
	errorStrategy 'retry'
	maxRetries 2
	// containerOptions '--memory=2g'
	cpus 1

	input:
	tuple val(timepoint), file(r1_fastq), file(r2_fastq), file(primers)

	output:
	tuple val(timepoint), file("*.sam"), file(primers)

	script:
	"""
	minimap2 -ax sr -t 1 -o ${timepoint}.sam ${params.refseq} ${r1_fastq} ${r2_fastq} && \
	rm -f `realpath *.fastq.gz`
	"""

}


process ILL_ALIGNMENT_POLISHING {

	// Using SAMTOOLS to clean and prep alignments

	tag "${timepoint}"
	
	publishDir params.results_data_files, pattern: '*.bam', mode: 'copy'
	publishDir params.results_data_files, pattern: '*.bam.bai', mode: 'move'

	input:
	tuple val(timepoint), file(sam), file(primers)

	output:
	tuple val(timepoint), file("*.bam"), emit: bam
	path("*bam.bai"), emit: bai

	script:
	"""

	cat ${sam} \
	| samtools ampliconclip -b ${primers} - \
	| samtools view -Sb - \
	| samtools sort - > ${timepoint}.bam && \
	rm -f `realpath ${sam}`

	samtools index ${timepoint}.bam

	"""

}

process ILL_LOWCOV_ANNOTATION {

	// Creating BED files that annotate regions of each genome with less than 20x coverage

	tag "${timepoint}"
	
	publishDir params.results_data_files, mode: 'copy'

	input:
	tuple val(timepoint), file(bam)

	output:
	file("*.bed")

	script:
	"""
	covtobed --max-cov=20 ${bam} > ${timepoint}.bed && \
	rm -f `realpath ${bam}`
	"""

}


process FIGURE_2A_PLOTTING {

	// Plotting consensus mutations

	publishDir params.visuals, pattern: '*.pdf', mode: 'move'
	publishDir params.results, pattern: '*.csv', mode: 'move'

	input:
	file(consensus)
	file(ont_bed_list)
	file(ill_bed_list)
	file(variant_table_list)

	output:
	val("finished"), emit: cue
	file("*.pdf")
	file("*.csv")

	script:
	"""
	Figure2A_consensus_mutations_through_time.R ${consensus} && \
	rm -f *.bed && \
	rm -f *.tsv
	"""

}


process TAR_FIG2A_DATA {
	
	publishDir params.results, pattern: '*.tar.xz', mode: 'move'
	
	cpus 4
	
	when:
	cue == "finished"
	
	input:
	val(cue)
	
	output:
	path("*.tar.xz")
	
	script:
	date = new java.util.Date().format('yyyyMMdd')
	
	"""
	tar -cf - ${params.results_data_files} | xz -9ve --threads=0 -z - > fig2a_data_${date}.tar.xz && \
	rm -rf ${params.results_data_files}
	"""
	
}


process FIGURE_2C_PLOTTING {

	// Plotting neutralization assay results

	publishDir params.visuals, pattern: '*.pdf', mode: 'move'

	output:
	file("*.pdf")

	script:
	"""
	Figure2C_neutralization_assay.R ${params.neut_data}
	"""

}

