#!/bin/bash -eux

conda activate variant-calling # see "SupplementalFigure1_conda_environment_setup.sh" for instructions on setting up an identical conda environment.

WD=${1:-/Volumes/working_ssd/e484t_manuscript/data/all_timepoints_raw_reads}
REF=${2:-/Volumes/working_ssd/e484t_manuscript/data/all_timepoints_raw_reads/NC_045512.2.fasta}


cd $WD

### TRIMMMING PAIRED READS ###
for i in $(ls *Illumina_R1*);
do
NAME=${i/_R1.fastq.gz/}
READS2=$(ls ${NAME}_R2*)
trimmomatic PE $i $READS2 \
${NAME}_R1_trimmed.fastq.gz ${NAME}_R1_unpaired.fastq.gz \
${NAME}_R2_trimmed.fastq.gz ${NAME}_R2_unpaired.fastq.gz \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done



### TRIMMING UNPAIRED READS ###
for i in $(ls *ONT*);
do
NAME=${i/.fastq.gz/}
trimmomatic SE $i ${NAME}_trimmed.fastq.gz \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 -phred64
done

for i in $(ls *IonTorrent*);
do
NAME=${i/.fastq.gz/}
trimmomatic SE $i ${NAME}_trimmed.fastq.gz \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done



### ALIGNING LONG READS ###
minimap2 -d ref.mmi $REF # indexing reference
for i in $(ls *ONT*);
do
NAME=${i/.fastq.gz/}
minimap2 -ax map-ont $REF $i | samtools view -Sb - | samtools sort - > ${NAME}.bam 
done



### ALIGNING SINGLE-END SHORT READS (ION TORRENT) ###
for i in $(ls *IonTorrent*);
do
NAME=${i/.fastq.gz/}
minimap2 -ax sr $REF $i | samtools view -Sb - | samtools sort - > ${NAME}.bam 
done



### ALIGNING PAIRED-END SHORT READS ###
for i in $(ls *Illumina_R1_trimmed*);
do
NAME=${i/_R1_trimmed.fastq.gz/}
READS2=$(ls ${NAME}_R2_trimmed*)
minimap2 -ax sr $REF $i $READS2 | samtools view -Sb - | samtools sort - > ${NAME}.bam 
done

for i in $(ls *bam);
do
samtools index $i
done



### VARIANT CALLING ###
# calling with bbmap callvariants.sh
bbmap=/Users/nicholasminor/bbmap
data=/Volumes/working_ssd/e484t_manuscript/data/all_timepoints_raw_reads
bamlist=$data/bam-list
DATE=$(date +'%Y%m%d')

cd $data
$bbmap/callvariants.sh list=$bamlist out=$data/alltimepoints_variants_${DATE}.vcf \
ref=$REF samstreamer=t ss=4 clearfilters multisample=t \
ploidy=1 mincov=0 overwrite=t

gunzip individual*.vcf.gz



### COLLECTING AVERAGE COVERAGE ###
for i in $(ls *bam);
do
bedtools genomecov -d -ibam $i | awk 'BEGIN {sum=0}; {sum+=$3}; END{print sum/NR}' < ${i}.avg_cov
done
