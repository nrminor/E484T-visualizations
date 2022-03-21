#!/bin/bash


PRIMERS=${1:-ref/primers_per_sample.txt}
DATE=$(date +'%Y%m%d')



## CLIPPING TO AMPLICONS AND CONVERTING TO BAMS
find "data/raw_reads" -maxdepth 1 -type f -name "*.sam" > data/raw_reads/sam_list.txt 

for i in `cat data/raw_reads/sam_list.txt `;

do
  f=$(basename "$i")
  NAME=${f/.sam/}
  primer_set=$(grep $i ref/primers_per_sample.txt | cut -f 2)
  primer_set=$(grep $NAME $PRIMERS | cut -f 2)
  cat $i \
  | samtools ampliconclip -b ref/$primer_set - \
  | samtools view -Sb - \
  | samtools sort - > data/raw_reads/${NAME}.bam
  samtools index data/raw_reads/${NAME}.bam
done
rm data/raw_reads/sam_list.txt
rm data/raw_reads/*.sam
