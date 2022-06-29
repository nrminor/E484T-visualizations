#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Define some helpful runtime parameters
params.primer_key = "$baseDir/config/fig2b_raw_read_guide.csv"
params.refdir = "$baseDir/resources"
params.refseq = "$baseDir/resources/reference.fasta"
params.refgff = "$baseDir/resources/MN9089473.gff3"
params.fasta_sep = "$baseDir/scripts/fasta_sep.R"
params.Ct_script = "$baseDir/scripts/Figure1A_Ct_through_time.R"
params.Ct_data = "$baseDir/data/Ct_timeline.csv"
params.cons_mutations = "$baseDir/scripts/Figure2A_consensus_mutations_through_time.R"
params.neut_script = "$baseDir/scripts/Figure2C_neutralization_assay.R"
params.neut_data = "$baseDir/data/antibody_potency.csv"
params.results = "$baseDir/results"
params.visuals = "$baseDir/results/visuals"


// Defining workflow for generating non-supplementary figures for Halfmann et al. 2022
workflow {

	// INPUT CHANNELS
	ch_consensus_seqs = Channel
		.fromPath('data/alltimepoints_20220222.fasta')

	ch_ont_reads = Channel // this channel streams in data on the timepoints sequenced on ONT instruments
		.fromPath ( "config/fig2b_raw_read_guide.csv" )
		.splitCsv ( sep: ",", header: true )
		.map      { row -> tuple(row.sra_id, row.file_basename, row.platform, file(row.primer_set)) }
		.filter   { it[2] == "ont" }

	ch_illumina_reads = Channel // this channel streams in data on the timepoints sequenced on Illumina instruments
		.fromPath ( "config/fig2b_raw_read_guide.csv" )
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
		pangolin --outfile 'prolonged_infection_lineage_report.csv' $consensus
		"""
}


process FIGURE_1A_PLOTTING {

	// Plotting Ct values

	publishDir params.visuals, pattern: '*.pdf', mode: 'move'

	output:
		file("*.pdf")

	script:
		"""
		Rscript ${params.Ct_script} ${params.Ct_data}
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
		Rscript ${params.fasta_sep} ${consensus}
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
	publishDir "results/data", mode: "copy"

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
	publishDir "data/reads", pattern: '*.fastq.gz'

	input:
		tuple val(sra_id), val(timepoint), val(platform), file(primers)

	output:
		tuple val(timepoint), file("*.fastq.gz"), file(primers)

	script:
		"""

		fasterq-dump ${sra_id} --concatenate-reads --skip-technical --quiet
		gzip ${sra_id}.fastq
		mv ${sra_id}.fastq.gz ${timepoint}.fastq.gz

		"""

}


process ONT_READ_MAPPING {

	// Mapping ONT reads from each timepoint

	tag "${timepoint}"
	publishDir "data/", pattern: '*.sam'

	input:
		tuple val(timepoint), file(fastq), file(primers)

	output:
		tuple val(timepoint), file("*.sam"), file(primers)

	script:
		"""
		minimap2 -ax map-ont ${params.refseq} ${fastq} > ${timepoint}.sam
		"""

}


process ONT_ALIGNMENT_POLISHING {

	// Using SAMTOOLS to clean and prep alignments

	tag "${timepoint}"
	publishDir "data/"

	input:
		tuple val(timepoint), file(sam), file(primers)

	output:
		tuple val(timepoint), file("*.bam"), emit: bam
		file("*.bam.bai")

	script:
		"""

		cat ${sam} \
		| samtools ampliconclip -b ${primers} - \
		| samtools view -Sb - \
		| samtools sort - > ${timepoint}.bam

		samtools index ${timepoint}.bam

		"""

}

process ONT_LOWCOV_ANNOTATION {

	// Creating BED files that annotate regions of each genome with less than 20x coverage

	tag "${timepoint}"
	publishDir "results/data", mode: "copy"

	input:
	tuple val(timepoint), file(bam)

	output:
	file("*.bed")

	script:
		"""
		covtobed --max-cov=20 ${bam} > ${timepoint}.bed
		"""

}


process GET_ILL_READS {

	tag "${timepoint}"
	publishDir "data/reads", pattern: '*.fastq.gz'

	input:
		tuple val(sra_id), val(timepoint), val(platform), file(primers)

	output:
		tuple val(timepoint), file("*_R1.fastq.gz"), file("*_R2.fastq.gz"), file(primers)

	script:
		"""

		fasterq-dump ${sra_id} --split-files --skip-technical --quiet
		gzip ${sra_id}_1.fastq ; mv ${sra_id}_1.fastq.gz ${timepoint}_R1.fastq.gz
		gzip ${sra_id}_2.fastq ; mv ${sra_id}_2.fastq.gz ${timepoint}_R2.fastq.gz

		"""

}


process ILL_READ_MAPPING {

	// Mapping ONT reads from each timepoint

	tag "${timepoint}"
	publishDir "data/", pattern: '*.sam'

	input:
		tuple val(timepoint), file(r1_fastq), file(r2_fastq), file(primers)

	output:
		tuple val(timepoint), file("*.sam"), file(primers)

	script:
		"""
		minimap2 -ax sr ${params.refseq} ${r1_fastq} ${r2_fastq} > ${timepoint}.sam
		"""

}


process ILL_ALIGNMENT_POLISHING {

	// Using SAMTOOLS to clean and prep alignments

	tag "${timepoint}"
	publishDir "data/"

	input:
		tuple val(timepoint), file(sam), file(primers)

	output:
		tuple val(timepoint), file("*.bam"), emit: bam
		file("*bam.bai")

	script:
		"""

		cat ${sam} \
		| samtools ampliconclip -b ${primers} - \
		| samtools view -Sb - \
		| samtools sort - > ${timepoint}.bam

		samtools index ${timepoint}.bam

		"""

}

process ILL_LOWCOV_ANNOTATION {

	// Creating BED files that annotate regions of each genome with less than 20x coverage

	tag "${timepoint}"
	publishDir "results/data", mode: "copy"

	input:
		tuple val(timepoint), file(bam)

	output:
		file("*.bed")

	script:
		"""
		covtobed --max-cov=20 ${bam} > ${timepoint}.bed
		"""

}


process FIGURE_2A_PLOTTING {

	// Plotting consensus mutations

	publishDir params.visuals, pattern: '*.pdf', mode: 'move'
	publishDir 'results/data/', pattern: '*.csv', mode: 'move'

	input:
		file(consensus)
		file(ont_bed_list)
		file(ill_bed_list)
		file(variant_table_list)

	output:
		file("*.pdf")
		file("*.csv")

	script:
		"""
		Rscript ${params.cons_mutations} ${consensus}
		"""

}


process FIGURE_2C_PLOTTING {

	// Plotting neutralization assay results

	publishDir params.visuals, pattern: '*.pdf', mode: 'move'

	output:
		file("*.pdf")

	script:
		"""
		Rscript ${params.neut_script} ${params.neut_data}
		"""

}

