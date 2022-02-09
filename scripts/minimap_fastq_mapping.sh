#!/bin/bash -eux



### SETUP ###
### ----- ###
REF=${1:-data/reference.fasta}
DATE=$(date +'%Y%m%d')



### ONT READS ###
### --------- ###
find "." -maxdepth 1 -type f -name "*ONT.fastq.gz" > ont_fastq_list.txt 

for i in `cat ont_fastq_list.txt `;
do
  f=$(basename "$i")
  NAME=${f/.fastq.gz/}
  minimap2 -ax map-ont $REF $f > ${NAME}.sam
done

rm ont_fastq_list.txt



### IONTORRENT READS ###
### ---------------- ###
find "." -maxdepth 1 -type f -name "*IonTorrent.fastq.gz" > ion_fastq_list.txt 

for i in `cat ion_fastq_list.txt `;
do
  f=$(basename "$i")
  NAME=${f/.fastq.gz/}
  minimap2 -ax sr $REF $f > ${NAME}.sam
done

rm ion_fastq_list.txt



### ILLUMINA READS ###
### ---------------- ###
find "." -maxdepth 1 -type f -name "*Illumina_R1.fastq.gz" > ill_fastq_list.txt 

for i in `cat ill_fastq_list.txt `;
do
  f=$(basename "$i")
  NAME=${f/_R1.fastq.gz/}
  READS2=${NAME}_R2.fastq.gz
  minimap2 -ax sr $REF $f $READS2 > ${NAME}.sam
done

rm ill_fastq_list.txt
