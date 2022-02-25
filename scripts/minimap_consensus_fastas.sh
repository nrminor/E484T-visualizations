#!/bin/bash


find "data/tmp" -maxdepth 1 -type f -name "*.fasta" > data/tmp/fasta_list.txt

for i in `cat data/tmp/fasta_list.txt  `;
do
  f=$(basename "$i")
  NAME=${f/.fasta/}
  minimap2 -a ref/reference.fasta $i > data/tmp/${NAME}.sam
done

rm data/tmp/fasta_list.txt
