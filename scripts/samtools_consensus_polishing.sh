#!/bin/bash

for i in `cat data/tmp/sam_list.txt `;
do
  f=$(basename "$i")
  NAME=${f/.sam/}
  cat $i \
  | samtools view -Sb - \
  | samtools sort - > data/tmp/${NAME}.bam
done

