#!/bin/bash -eux

DATE=$(date +'%Y%m%d')

find "." -maxdepth 1 -type f -name "*.mpileup" > mpileup_list.txt 

for i in `cat mpileup_list.txt `;
do
  f=$(basename "$i")
  NAME=${f/.mpileup/}
  cat $f \
  | ivar consensus -p ${NAME}_ivar_${DATE} -t 0.5 -m 20
done

rm mpileup_list.txt
rm *.mpileup

