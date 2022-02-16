#!/bin/bash


find "data/b12_enriched_global/tmp" -maxdepth 1 -type f -name "*.fasta" > data/b12_enriched_global/tmp/fasta_list.txt

for i in `cat data/b12_enriched_global/tmp/fasta_list.txt  `;
do
  f=$(basename "$i")
  NAME=${f/.fasta/}
  minimap2 -a ref/reference.fasta $i > data/b12_enriched_global/tmp/${NAME}.sam
done

find "data/b12_enriched_global/tmp" -maxdepth 1 -type f -name "*.sam" > data/b12_enriched_global/tmp/sam_list.txt
