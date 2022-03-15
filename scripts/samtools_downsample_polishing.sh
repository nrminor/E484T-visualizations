#!/bin/bash


PRIMERS=${1:-ref/primers_per_sample.txt}
DATE=$(date +'%Y%m%d')



# RE-SORTING & CLEANING UP DOWNSAMPLED READS
for i in `cat data/raw_reads/bam_list.txt `;
do
  f=$(basename "$i")
  NAME=${f/.bam/}
  amplicons=$(grep $NAME $PRIMERS | cut -f 3)
  samtools sort data/raw_reads/${NAME}_downsampled.bam > data/raw_reads/${NAME}.bam
  rm data/raw_reads/${NAME}_downsampled.bam ; rm data/raw_reads/${NAME}_collated.bam
  samtools index data/raw_reads/${NAME}.bam
done

rm data/raw_reads/bam_list.txt
