#!/bin/bash

DATE=$(date +'%Y%m%d')

find "data/tmp" -maxdepth 1 -type f -name "*.mpileup" > data/tmp/mpileup_list.txt 

for i in `cat data/tmp/mpileup_list.txt `;
do
  f=$(basename "$i")
  NAME=${f/.mpileup/}
  cat $i \
  | ivar variants -p data/tmp/${NAME}_consensus_variant_table -t 0 -m 1 -q 1 -r ref/reference.fasta -g ref/MN9089473.gff3
done

rm data/tmp/mpileup_list.txt
rm data/tmp/*.mpileup

