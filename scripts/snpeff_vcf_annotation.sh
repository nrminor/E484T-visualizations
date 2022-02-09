#!/bin/bash -eux

WD=$(pwd)

find "." -maxdepth 1 -type f -name "*.vcf" > vcf_list.txt

for i in `cat vcf_list.txt`;
do
  # get the basename as i used full path in find (which I typically find more usefull, either way a ./ will be added as a prefix with find.)
  f=$(basename "$i")
  NAME=${f/.vcf/}
  cat $f \
  | snpEff -c $WD/snpEff_SARS-CoV-2_database/MN908947.3.config -onlyProtein MN908947.3 > ${NAME}.ann.vcf
done

rm vcf_list.txt
