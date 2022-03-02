#!/bin/bash

WD=${1:-data/raw_reads} # relative path to FASTAs

find $WD -maxdepth 1 -type f -name "*.fasta" > $WD/fasta_list.txt

for i in `cat $WD/fasta_list.txt `;
do
  f=$(basename "$i")
  NAME=${f/.fasta/}
  minimap2 -a ref/reference.fasta $i > $WD/${NAME}.sam
done

rm $WD/fasta_list.txt
