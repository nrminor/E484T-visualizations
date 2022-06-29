## Introduction

This workflow was created to efficiently gather, process, and visualize data presented in [Halfmann et al. 2022, _Evolution of a globally unique SARS-CoV-2 Spike E484T monoclonal antibody escape mutation in a persistently infected, immunocompromised individual_](https://www.medrxiv.org/content/10.1101/2022.04.11.22272784v1). Raw reads from the course of the prolonged infection are pulled automatically from SRA BioProject PRJNA836936. Consensus sequences are bundled together with the workflow documents and are included in [the project GitHub repository](https://github.com/dholab/E484T-visualizations/tree/main). Also included is a table of neutralization assay results and a table of qPCR cycle threshold (Ct) values from the infection.

Note that this workflow does not include processing behind Halfmann et al. 2022 supplemental tables and figures. These workflows are bundled separately and can be accessed by following the links in these GitHub repos:
- [Supplemental Figure 1](https://github.com/nrminor/prolonged-infection-suppfig1) - (still under development)
- [Supplemental Figure 2](https://github.com/nrminor/prolonged-infection-suppfig2)

## Getting started

To reproduce our results, we recommend users simply clone the repository to their machine of choice with `git clone https://github.com/dholab/E484T-visualizations.git`. Once the repository is in place, users should make sure the [NextFlow](https://www.nextflow.io/) workflow manager is installed. If it isn't installed, go through the following steps:

1. Install the miniconda python distribution (https://docs.conda.io/en/latest/miniconda.html)
2. Install the `mamba` package installation tool:
   `conda install -y -c conda-forge mamba`
3. Install `nextflow` and other packages needed for the workflow by changing into the workflow directory and running:
   ` mamba env create -f config/envs/prolonged_infection.yaml -n prolonged-infection `
4. Activate the new environment:
   ` conda activate prolonged-infection `

To double check that the installation was successful, type `nextflow -v` into the terminal. If it returns something like `nextflow version 21.04.0.5552`, you are set and ready to proceed.

To run the workflow, simply change into the workflow directory and run the following in the BASH terminal:

```
nextflow run prolonged_infection_workflow.nf
```

If the workflow runs partway, but a computer outage or other issue interrupts its progress, no need to start over! Instead, run:

```
nextflow run prolonged_infection_workflow.nf -resume
```

Finally, if you're like me and you appreciate additional information about how the workflow ran, run the workflow with any or all of the following flags:

```
nextflow prolonged_infection_workflow.nf -with-report -with-timeline -with-dag prolonged_infection_dag.png
```

Note that NextFlow's DAG-plotting requires that you have install GraphViz, which is easiest to do via the intructions on [GraphViz's website](https://graphviz.org/download/).

### Configuration

The following runtime parameters have been set for the whole workflow:

- `primer_key` - path to an important CSV that lists which primer set was used for each sample
- `refdir` - path to the directory where the reference sequence and primer BED files are located
- `refseq` - path to the SARS-CoV-2 Wuhan-1 sequence from GenBank Accession MN9089473.3
- `refgff` - path to the gene feature file for the SARS-CoV-2 Wuhan-1 sequence
- `fasta_sep` - path to a script that separates out our consensus sequences into individual FASTA files. This script was necessary because GISAID consensus sequence strain names, which I use to name each FASTA, contain backslashes. The script - `fasta_sep.R` parses the string for the GISAID strain name and replaces slashes with underscores before naming the FASTA.
- `Ct_script` - path to the script that plots qPCR Ct values through time
- `Ct_data` - path to the CSV file of Ct values
- `cons_mutations` - path to the script that plots Figure 2A
- `neut_script` - path to the script that plots neutralization assay results
- `neut_data` - path to the CSV file of neutralization assay data
- `results` - path to the results directory, the workflow's default output location
- `results_data_files` - path to a subdirectory of `results/` where data files like VCFs are placed.
- `visuals` - path to a subdirectory of `results/` where graphics are stored. This is where the final figure PDFs are placed.

These parameters can be altered in the command line with a double-dash flag, like so:

```
nextflow run prolonged_infection_workflow.nf --results ~/Downloads
```

### Bundled data

Bundled together with the workflow are:

- primer BED files for ARTIC v3, ARTIC v4.1, and MIDNIGHT v1 are in `ref/`
- the SARS-CoV-2 Wuhan-1 sequence from GenBank Accession MN9089473.3 is in `ref/` and is called reference.fasta
- the MN9089473.3 gene coordinate/codon GFF file is also in `ref/` and is called MN9089473.gff3
- a CSV file of qPCR cycle threshold values through the course of the prolonged infection is in `data/`
- a CSV file of neutralization assay results is also in `data/`
- a very important file called `fig2b_raw_read_guide.csv` is in `config/`. This file enables the workflow to use the correct primers for each sample in `samtools ampliconclip`.

## Workflow summary

- the first process identifies the SARS-CoV-2 pango lineage for the consensus sequences from each timepoint
- qPCR cycle threshold values are then plotted
- each sequence from the FASTA of consensus sequences is isolated into its own FASTA file
- the consensus FASTAs are converted to SAM format (mapped to the Wuhan-1 SARS-CoV-2 sequence)
- a samtools `mpileup` is generated for each consensus SAM
- the consensus pileups are piped into iVar, which variant calls each consensus sequence and also identifies protein effects
- SRA accession numbers are taken from `config/fig2b_raw_read_guide.csv` and used to pull FASTQs for each sequencing timepoint. Because NCBI's `fastq-dump` is very slow, this step takes the longest. _NOTE_ that because the settings for mapping these reads depend on sequencing platform (e.g., Illumina vs. Oxford Nanopore), there are near-identical processes laid out for reads from each sequencing platform.
- the reads are mapped to the Wuhan-1 SARS-CoV-2 sequence
- the mapped reads in SAM format are then clipped down to amplicons and converted to BAM format
- a tool called `covtobed` then identifies regions in any of these BAMs where depth of coverage is less than 20 reads, indicating an amplicon dropout
- the amplicon dropout BED files, the consensus sequence FASTA, and the iVar variant tables are then plotted with an R script
- the neutralization assay results are plotted with an R script

## Output files

- `data/reads/[timepoint_day_sequencingplatform].fastq.gz` - Raw reads for each timepoint that we have deposited in Sequence Read Archive (SRA) BioProject PRJNA836936
- `data/[timepoint_day_sequencingplatform].sam` - reads aligned to the SARS-CoV-2 Wuhan-1 sequence in human-readable SAM format
- `data/[timepoint_day_sequencingplatform].bam` - reads aligned to the SARS-CoV-2 Wuhan-1 sequence, filtered down to only reads that are within the span of each amplicon in the timepoint's primer set, and converted to the smaller, non-human readable BAM format.
- `results/data/[timepoint_day_sequencingplatform].bed` - BED files showing genome regions where fewer than 20 reads aligned. These low-depth regions are likely to be regions where amplification failed, i.e., amplicon dropouts.
- `results/data/[GISAID_strain_name]_consensus_variant_table.tsv` - iVar tables of single-nucleotide variants identified in the consensus sequence of each timepoint, along with codon information and protein effects for each variant.
- `results/data/patient_variant_counts.csv` - simple table recording the number of variants iVar identified at each timepoint. This is used as an input for the script that generates Supplementary Figure 2.
- `results/prolonged_infection_lineage_report.csv` - pango lineages classified for each timepoint.
- `results/visuals/fig1a_ct_values_thru_time.pdf` - PDF file of Figure 1A, which can be polished in vector graphics software like Adobe Illustrator.
- `results/visuals/fig2a_consensus_mutations.pdf` - PDF file of Figure 2A, which can be polished in vector graphics software like Adobe Illustrator.
- `results/visuals/fig2c_neutralization_assay.pdf` - PDF file of Figure 2C, which can be polished in vector graphics software like Adobe Illustrator.

Note that figures 1B and 2B are produced manually.

## For more information

This workflow was created by Nicholas R. Minor. To report any issues, please visit the [GitHub repository for this project](https://github.com/dholab/E484T-visualizations/tree/main).

### NextFlow Learning Resources

- [NBI Sweden Reproducibility Workshop GitHub](https://github.com/NBISweden/workshop-reproducible-research/tree/main/pages/nextflow) - This one is more up to date than the version on canva, which is the next link
- [NBI Sweden Reproducibility Workshop Canva](https://uppsala.instructure.com/courses/51980/pages/nextflow-1-introduction?module_item_id=328997)
- [NextFlow Official Documentation](https://www.nextflow.io/docs/latest/index.html)
- [Exhaustive NextFlow Github.io Tutorial](https://nextflow-io.github.io/patterns/index.html#_basic_patterns)
- [Seqera Labs NextFlow Tutorial](https://training.seqera.io/)
- [Tronflow NextFlow tutorial](https://tronflow-docs.readthedocs.io/en/latest/02_tutorial.html)
- [Google Cloud NextFlow tutorial](https://cloud.google.com/life-sciences/docs/tutorials/nextflow)
- [Data Carpentries Incubator Workshop on NextFlow](https://carpentries-incubator.github.io/workflows-nextflow/index.html) - This one is especially thorough and helpful, and also introduces nf-core.
- [NextFlow Cheatsheet PDF](https://github.com/danrlu/Nextflow_cheatsheet/blob/main/nextflow_cheatsheet.pdf) - Not a tutorial, but this is a useful PDF to have on hand.
