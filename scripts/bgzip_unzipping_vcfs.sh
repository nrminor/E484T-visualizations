find "." -maxdepth 1 -type f -name "individual*.vcf.gz" > vcf_list.txt 

for i in `cat vcf_list.txt `;
do
  f=$(basename "$i")
  NAME=${f/.vcf.gz/}
  bgzip -d -f $f
done
