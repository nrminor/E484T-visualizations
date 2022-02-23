#!/bin/bash

find "data/" -maxdepth 1 -type f -name "*.vcf" > data/vcf_list.txt

for i in `cat data/vcf_list.txt`;
do
  # get the basename as i used full path in find (which I typically find more usefull, either way a ./ will be added as a prefix with find.)
  f=$(basename "$i")
  NAME=${f/.vcf/}
  cat $i \
  | snpEff -c data/snpEff_SARS-CoV-2_database/MN908947.3.config -onlyProtein MN908947.3 > data/${NAME}.ann.vcf
done

rm data/vcf_list.txt
