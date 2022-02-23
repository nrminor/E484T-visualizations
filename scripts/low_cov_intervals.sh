#!/bin/bash



### SETUP ###
### ----- ###
READS=${1:-data/raw_reads} # relative path to BAMs
OUTPUT=${2:-data/low_coverage_regions} # relative path to desired output directory

find $READS -maxdepth 1 -type f -name "*ONT.bam" > $READS/ont_bam_list.txt 

for i in `cat $READS/ont_bam_list.txt  `;
do
  f=$(basename "$i")
  NAME=${f/.bam}
  covtobed --max-cov=20 $i > $OUTPUT/${NAME}.bed
done

rm $READS/ont_bam_list.txt
