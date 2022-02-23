find "data/" -maxdepth 1 -type f -name "individual*.vcf.gz" > data/vcf_list.txt 

for i in `cat data/vcf_list.txt `;
do
  f=$(basename "$i")
  NAME=${f/.vcf.gz/}
  bgzip -d -f $i
done

rm data/vcf_list.txt
