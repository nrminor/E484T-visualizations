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


#### FIGURE 1A ####
#### --------- ####
cd $workingdir
Rscript $scripts/Figure1A_Ct_through_time.R $workingdir



#### FIGURE 2B ####
#### --------- ####
# will add bedtools code to replace Geneious low-coverage region finding
Rscript $scripts/Figure2A_consensus_mutations_through_time.R $workingdir



#### FIGURE 2C ####
#### --------- ####
Rscript $scripts/Figure2C_neutralization_assay.R $workingdir



#### SUPPLEMENTAL FIGURE 1 ####
#### --------------------- ####

#### MAPPING READS ####
# Indexing the reference fasta, mapping reads, cutting out primers, converting, and sorting alignments for each sample
READS=$data/raw_reads
REF=$workingdir/ref/reference.fasta
PRIMERS=$workingdir/ref/ARTICv3_primers.bed

cp $REF $READS
cp $PRIMERS $READS

cd $READS
REF=$(basename "$REF")
PRIMERS=$(basename "$PRIMERS")
DOCKER_RUN=(docker run -it --user $(id -u):$(id -g) -v $(pwd)/:/scratch -w /scratch)

cp $scripts/minimap_fastq_mapping.sh .
$DOCKER_RUN biocontainers/minimap2:v2.15dfsg-1-deb_cv1 \
/bin/bash minimap_fastq_mapping.sh $REF
rm minimap_fastq_mapping.sh

cp $scripts/samtools_alignment_polishing.sh .
$DOCKER_RUN quay.io/biocontainers/samtools:1.14--hb421002_0 \
/bin/bash samtools_alignment_polishing.sh $PRIMERS
rm samtools_alignment_polishing.sh

#### CONSENSUS- AND VARIANT-CALLING ####
# Calling a consensus for each sample with iVar
cp $scripts/ivar_consensus_calling.sh .
$DOCKER_RUN andersenlabapps/ivar:1.3.1 \
/bin/bash ivar_consensus_calling.sh
rm ivar_consensus_calling.sh

# Calling variants with BBTools
find "." -maxdepth 1 -type f -name "*.bam" > bam_list.txt 
$DOCKER_RUN quay.io/biocontainers/bbmap:38.93--he522d1c_0 \
callvariants.sh list=bam_list.txt out=alltimepoints_variants_${DATE}.vcf \
ref=$REF samstreamer=t ss=4 multisample=t clearfilters \
ploidy=1 mincov=0 minallelefraction=0.002 overwrite=t
rm bam_list.txt 

# Pangolin lineage identification for our consensus # 3.1.19-pangolearn-2022-01-20
cat *${DATE}.fa > ivar_consensus_seqs_${DATE}.fasta
rm *${DATE}.fa
$DOCKER_RUN --pull always staphb/pangolin:latest \
pangolin --outfile ivar_lineage_report.csv ivar_consensus_seqs_${DATE}.fasta

# tidying up
rm $REF
rm $PRIMERS
mv individual* ..
cp ivar_consensus_seqs_${DATE}.fasta ..
cp ivar_lineage_report.csv ..


#### ANNOTATING BBTOOLS VCFs ####
cd ..
DOCKER_RUN=(docker run -it --user $(id -u):$(id -g) -v $(pwd)/:/scratch -w /scratch)
mv ivar_consensus_seqs_${DATE}.fasta alltimepoints_${DATE}.fasta
cp $scripts/bgzip_unzipping_vcfs.sh .
docker run -it --user $(id -u):$(id -g) -v $(pwd)/:/scratch -w /scratch bioslimcontainers/tabix:1.7 \
/bin/bash bgzip_unzipping_vcfs.sh
rm bgzip_unzipping_vcfs.sh

cp $scripts/snpeff_vcf_annotation.sh .
docker run -it --user $(id -u):$(id -g) -v $(pwd)/:/data -w /data bioinfoipec/snpeff:latest \
/bin/bash snpeff_vcf_annotation.sh
rm snpeff_vcf_annotation.sh
rm *ONT.vcf
rm *Illumina.vcf
rm *IonTorrent.vcf

# Rscript $scripts/SupplementalFigure1_allele_frequency_plot.R $workingdir # HIGHLY RECOMMENDED TO RUN THIS SEPARATELY; IT WILL TAKE A DAY OR TWO



#### SUPPLEMENTAL FIGURE 2 ####
#### --------------------- ####
# git clone https://github.com/nextstrain/ncov.git $scripts/ncov
bash $scripts/SupplementalFigure2_gisaid_reformatting_subsampling.sh \
	-d $data/b12_enriched_global -s $data/b12_enriched_global/sequences_fasta_2022_02_15.tar.xz \
	-m $data/b12_enriched_global/metadata_tsv_2022_02_15.tar.xz \
	-i $data/b12_enriched_global/UShER_b12_relatives_metadata.tsv \
	-p b12_enriched_global
Rscript $scripts/SupplementalFigure2_VOC_plot.R $workingdir
