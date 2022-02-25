#!/bin/bash

find "data/tmp" -maxdepth 1 -type f -name "*.sam" > data/tmp/sam_list.txt

for i in `cat data/tmp/sam_list.txt `;
do
  f=$(basename "$i")
  NAME=${f/.sam/}
  cat $i \
  | samtools view -Sb - \
  | samtools sort - > data/tmp/${NAME}.bam
  samtools mpileup -aa -A -d 600000 -B -Q 0 -f ref/reference.fasta --output data/tmp/${NAME}.mpileup data/tmp/${NAME}.bam
done

rm data/tmp/sam_list.txt
