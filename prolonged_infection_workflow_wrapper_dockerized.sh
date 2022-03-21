#!/bin/bash



#### SCRIPT DESCRIPTION ####
#### ------------------ ####
# This wrapper script will establish the file architecture for the prolonged infection manuscript, 
# and connect the necessary files to the scripts used to produce each of the analyses/visualizations 
# in the manuscript. It will pull all necessary Docker containers and take advantage of the scripts
# present in the project's data architecture. Eventually, this rudimentary shell wrapper will be replaced
# with a Snakemake workflow, which will make running the whole workflow more efficient and reproducible.
# NOTE: 
# This script is meant to be run on a linux-based system with R and Docker installed. Your machine
# must also have an up to date local file-tree database, which you can do with the following:
# For Macs, run:
# sudo launchctl load -w /System/Library/LaunchDaemons/com.apple.locate.plist
# sudo /usr/libexec/locate.updatedb
# On a Linux system, run:
# updatedb



#### DATA ARCHITECTURE ####
#### ----------------- ####
while getopts :w:d:s:v: flag
do
	case "${flag}" in
		w) workingdir=${OPTARG};;
		:) workingdir=$(locate prolonged_infection_workflow_wrapper_dockerized.sh | grep $(whoami) | sed 's/prolonged_infection_workflow_wrapper_dockerized.sh//') ;;
		d) data=${OPTARG};;
		:) data=$workingdir/data ;;
		s) scripts=${OPTARG};;
		:) scripts=$workingdir/scripts ;;
		v) visuals=${OPTARG};;
		:) visuals=$workingdir/visuals ;;
	esac
done

DATE=$(date +'%Y%m%d')
cd $workingdir
DOCKER_RUN=(docker run -it --user $(id -u):$(id -g) -v $(pwd)/:/scratch -w /scratch)



# Pangolin lineage identification for our consensus # 3.1.19-pangolearn-2022-01-20
$DOCKER_RUN --pull always staphb/pangolin:latest \
pangolin --outfile readables/ivar_lineage_report.csv data/alltimepoints_20220222.fasta



#### FIGURE 1A ####
#### --------- ####
Rscript $scripts/Figure1A_Ct_through_time.R $workingdir



#### FIGURE 2B ####
#### --------- ####

### CALLING VARIANTS FROM CONSENSUS SEQUENCES ###
mkdir $data/tmp/
# separating each consensus sequence into its own fasta
Rscript $scripts/fasta_sep.R $data $data/alltimepoints_20220222.fasta

# "converting" FASTAs into SAMs, then BAMS, and finally into pileup format for iVar
$DOCKER_RUN quay.io/biocontainers/minimap2:2.24--h5bf99c6_0 \
/bin/bash scripts/minimap_consensus_fastas.sh data/tmp
$DOCKER_RUN quay.io/biocontainers/samtools:1.14--hb421002_0 \
/bin/bash scripts/samtools_consensus_polishing.sh

# creating non-VCF, simplified variant tables with iVar
$DOCKER_RUN andersenlabapps/ivar:1.3.1 \
/bin/bash scripts/ivar_variant_tables.sh

# moving VCF into position for plotting and clearing tmp files
mv data/tmp/*variant_table.tsv data/
rm -r $data/tmp/

### DEFINING LOW-COVERAGE REGIONS FROM READS ###
# Mapping reads, cutting out primers, converting, downsampling, and sorting alignments for each sample
$DOCKER_RUN biocontainers/minimap2:v2.15dfsg-1-deb_cv1 \
/bin/bash scripts/minimap_fastq_mapping.sh
$DOCKER_RUN quay.io/biocontainers/samtools:1.14--hb421002_0 \
/bin/bash scripts/samtools_alignment_polishing.sh ref/primers_per_sample.txt

# defining low-coverage intervals in ONT data with "covtobed"
$DOCKER_RUN quay.io/biocontainers/covtobed:1.3.5--h36a6f06_0 \
/bin/bash scripts/low_cov_intervals.sh 

### PLOTTING THE DATA ###
Rscript $scripts/Figure2A_consensus_mutations_through_time.R $workingdir



#### FIGURE 2C ####
#### --------- ####
Rscript $scripts/Figure2C_neutralization_assay.R $workingdir



#### SUPPLEMENTAL FIGURE 1 ####
#### --------------------- ####

### VARIANT-CALLING MINOR READ VARIANTS ###
# Calling variants with BBTools
find "data/raw_reads/" -maxdepth 1 -type f -name "*Illumina.bam" > data/raw_reads/Illumina_bam_list.txt 
$DOCKER_RUN quay.io/biocontainers/bbmap:38.93--he522d1c_0 \
callvariants.sh list=data/raw_reads/Illumina_bam_list.txt out=data/alltimepoints_minor_variants_${DATE}.vcf \
ref=ref/reference.fasta samstreamer=t ss=4 multisample=t clearfilters \
ploidy=1 mincov=20 minallelefraction=0.002 callsub=t calldel=f callins=f nopassdot=t overwrite=t
rm data/raw_reads/Illumina_bam_list.txt 
mv individual*.vcf.gz $data
gunzip $data/individual*.vcf.gz

### ANNOTATING BBTOOLS VCFs ###
docker run -it --user $(id -u):$(id -g) -v $(pwd)/:/data -w /data bioinfoipec/snpeff:latest \
/bin/bash scripts/snpeff_vcf_annotation.sh
rm snpEff_*

Rscript $scripts/SupplementalFigure1_iSNV_frequency_plot.R $workingdir



#### SUPPLEMENTAL FIGURE 2 ####
#### --------------------- ####
# git clone https://github.com/nextstrain/ncov.git $scripts/ncov
ls -t -1 data/b12_enriched_global/sequences*.tar.xz > data/b12_enriched_global/sequence_tars.txt
sequence=$(basename $(head -n 1 data/b12_enriched_global/sequence_tars.txt)) && rm data/b12_enriched_global/sequence_tars.txt
ls -t -1 data/b12_enriched_global/metadata*.tar.xz > data/b12_enriched_global/metadata_tars.txt
metadata=$(basename $(head -n 1 data/b12_enriched_global/metadata_tars.txt)) && rm data/b12_enriched_global/metadata_tars.txt

bash $scripts/SupplementalFigure2_gisaid_reformatting_subsampling.sh \
	-d $workingdir/data/b12_enriched_global -s data/b12_enriched_global/$sequence \ # ALL PATHS MUST BE RELATIVE
	-m data/b12_enriched_global/$metadata \
	-i data/b12_enriched_global/UShER_b12_relatives_metadata.tsv \
	-p b12_enriched_global
Rscript $scripts/SupplementalFigure2_gisaid_roottotip_plot.R $workingdir
