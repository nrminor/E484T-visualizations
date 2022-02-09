#!/bin/bash -eux


PRIMERS=${1:-data/ARTICv3_primers.bed}
DATE=$(date +'%Y%m%d')

find "." -maxdepth 1 -type f -name "*.sam" > sam_list.txt 

for i in `cat sam_list.txt `;
do
  f=$(basename "$i")
  NAME=${f/.sam/}
  cat $f \
  | samtools ampliconclip -b $PRIMERS - \
  | samtools view -Sb - \
  | samtools sort - > ${NAME}.bam
  samtools index ${NAME}.bam
  samtools mpileup -aa -A -d 600000 -B -Q 0 --output ${NAME}.mpileup ${NAME}.bam
done

rm sam_list.txt
rm *.sam

