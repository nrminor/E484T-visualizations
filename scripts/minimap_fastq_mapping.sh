#!/bin/bash



### SETUP ###
### ----- ###
READS=${1:-data/raw_reads} # relative path to raw read FASTQs
REF=${2:-ref/reference.fasta}
DATE=$(date +'%Y%m%d')



### ONT READS ###
### --------- ###
find $READS -maxdepth 1 -type f -name "*ONT.fastq.gz" > $READS/ont_fastq_list.txt 

for i in `cat $READS/ont_fastq_list.txt `;
do
  f=$(basename "$i")
  NAME=${f/.fastq.gz/}
  minimap2 -ax map-ont $REF $i > $READS/${NAME}.sam
done

rm $READS/ont_fastq_list.txt



### ILLUMINA READS ###
### ---------------- ###
find $READS -maxdepth 1 -type f -name "*R1_Illumina.fastq.gz" > $READS/ill_fastq_list.txt 

for i in `cat $READS/ill_fastq_list.txt `;
do
  f=$(basename "$i")
  NAME=${f/_R1_Illumina.fastq.gz/}
  READS2=${NAME}_R2_Illumina.fastq.gz
  minimap2 -ax sr $REF $i $READS/$READS2 > $READS/${NAME}_Illumina.sam
done

rm $READS/ill_fastq_list.txt
