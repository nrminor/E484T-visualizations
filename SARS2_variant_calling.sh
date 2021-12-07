#!/bin/bash -eux

conda activate variant-calling

WD=${1:-/Volumes/working_ssd/e484t_manuscript/data/all_timepoints_raw_reads}
REF=${2:-/Volumes/working_ssd/e484t_manuscript/data/all_timepoints_raw_reads/NC_045512.2.fasta}


cd $WD

### TRIMMMING PAIRED READS ###
for i in $(ls *Illumina_R1*);
do
NAME=${i/_R1.fastq.gz/}
READS2=$(ls ${NAME}_R2*)
trimmomatic PE $i $READS2 \
${NAME}_R1_trimmed.fastq.gz ${NAME}_R1_unpaired.fastq.gz \
${NAME}_R2_trimmed.fastq.gz ${NAME}_R2_unpaired.fastq.gz \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done



### TRIMMING UNPAIRED READS ###
for i in $(ls *ONT*);
do
NAME=${i/.fastq.gz/}
trimmomatic SE $i ${NAME}_trimmed.fastq.gz \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 -phred64
done

for i in $(ls *IonTorrent*);
do
NAME=${i/.fastq.gz/}
trimmomatic SE $i ${NAME}_trimmed.fastq.gz \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done



### ALIGNING LONG READS ###
minimap2 -d ref.mmi $REF # indexing reference
for i in $(ls *ONT*);
do
NAME=${i/.fastq.gz/}
minimap2 -ax map-ont $REF $i | samtools view -Sb - | samtools sort - > ${NAME}.bam 
done



### ALIGNING SINGLE-END SHORT READS (ION TORRENT) ###
for i in $(ls *IonTorrent*);
do
NAME=${i/.fastq.gz/}
minimap2 -ax sr $REF $i | samtools view -Sb - | samtools sort - > ${NAME}.bam 
done



### ALIGNING PAIRED-END SHORT READS ###
for i in $(ls *Illumina_R1_trimmed*);
do
NAME=${i/_R1_trimmed.fastq.gz/}
READS2=$(ls ${NAME}_R2_trimmed*)
minimap2 -ax sr $REF $i $READS2 | samtools view -Sb - | samtools sort - > ${NAME}.bam 
done

for i in $(ls *bam);
do
samtools index $i
done


### VARIANT CALLING ###
# calling with bcftools
# DATE=$(date +'%Y%m%d')
# bcftools mpileup -a AD,DP,SP --threads 4 --max-idepth 1000000 --max-depth 1000000 \
# --fasta-ref ${REF} $(ls *bam) | bcftools call --threads 4 -o alltimepoints_variants_${DATE}.vcf) \
# -Ov --ploidy-file <(echo '* * * * 1') --keep-alts --variants-only --multiallelic-caller
# for i in $(bcftools query -l alltimepoints_variants_${DATE}.vcf);
# do
# bcftools view -s $i bcftools view -s > variants_${i}.vcf
# done
# 
# for i in $(ls *bam);
# do
# NAME=${i/.bam/}
# bcftools mpileup -a AD,DP,SP --threads 4 --max-idepth 1000000 --max-depth 1000000 \
# --fasta-ref $REF $i -Ou | bcftools call --threads 4 -o ${NAME}_variants_${DATE}.vcf \
# -Ov --ploidy-file <(echo '* * * * 1') --keep-alts --variants-only --multiallelic-caller
# done

# calling with Freebayes
# DATE=$(date +'%Y%m%d')
# ls *bam > bam-list
# freebayes-parallel <(fasta_generate_regions.py NC_045512.2.fasta.fai 1000) 4 \
# -f $REF --gvcf --ploidy 1 --max-coverage 1000000 \
# --bam-list bam-list > alltimepoints_variants_${DATE}_freebayes.vcf
# 
# for i in $(ls *bam);
# do
# NAME=${i/.bam/}
# freebayes-parallel <(fasta_generate_regions.py NC_045512.2.fasta.fai 1000) 4 \
# -f $REF --haplotype-length 0 --min-alternate-count 1 --ploidy 1 --max-coverage 1000000 \ 
# --min-alternate-fraction 0 --pooled-continuous --report-monomorphic \
# --gvcf -b $i > ${NAME}_variants_${DATE}_freebayes.vcf
# done
# 
# for i in $(ls *_freebayes.vcf);
# do
# NAME=${i/_freebayes.vcf/}
# sed -i '' -e "s/unknown/$NAME/g" "$i" # replace -i '' -e flags with -i alone when not running on a mac # https://stackoverflow.com/questions/19456518/error-when-using-sed-with-find-command-on-os-x-invalid-command-code
# bgzip -c $i > ${NAME}_freebayes.vcf.gz
# bcftools index ${NAME}_freebayes.vcf.gz
# done
# 
# # ls *freebayes.vcf.gz > vcf-list
# # bcftools merge --file-list vcf-list --gvcf $REF --merge both --threads 4 -O v -o alltimepoints_variants_${DATE}_freebayes.vcf
# 
# 
# # calling with VarDict # https://github.com/AstraZeneca-NGS/VarDict
# AF_THR="0.0001" # minimum allele frequency
# DATE=$(date +'%Y%m%d')
# for i in $(ls *bam);
# do
# NAME=${i/.bam/}
# vardict -G $REF -f $AF_THR -V $AF_THR -b $i -v ${NAME}_variants_${DATE}_vardict.vcf
# # ~or~ #
# # vardict -G $REF -f $AF_THR -V $AF_THR -b $i | var2vcf_valid.pl -N sample_name -E -f $AF_THR
# done
# 
# # calling with LoFreq
# ls *bam > bam-list
# samtools merge alltimepoints_merged.bam *.bam -o
# lofreq call-parallel --pp-threads 4 -f $REF -o ${NAME}_variants_${DATE}_lowfreq.vcf alltimepoints_merged.bam

# calling with bbmap callvariants.sh
bbmap=/Users/nicholasminor/bbmap
data=/Volumes/working_ssd/e484t_manuscript/data/all_timepoints_raw_reads
bamlist=$data/bam-list
DATE=$(date +'%Y%m%d')

cd $data
$bbmap/callvariants.sh list=$bamlist out=$data/alltimepoints_variants_${DATE}.vcf \
ref=$REF samstreamer=t ss=4 clearfilters multisample=t \
ploidy=1 mincov=0 overwrite=t

gunzip individual*.vcf.gz


### COLLECTING AVERAGE COVERAGE ###
for i in $(ls *bam);
do
bedtools genomecov -d -ibam $i | awk 'BEGIN {sum=0}; {sum+=$3}; END{print sum/NR}' < ${i}.avg_cov
done
