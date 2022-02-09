#!/bin/bash -eux

conda activate variant-calling # see "SupplementalFigure1_conda_environment_setup.sh" for instructions on setting up an identical conda environment.

WD=${1:-data/raw_reads}
REF=${2:-../reference.fasta}
bbmap=${3:-../../software/bbmap}


cd $WD

### TRIMMMING PAIRED READS ###
for i in $(ls *Illumina_R1.fastq.gz);
do
NAME=${i/_R1.fastq.gz/}
READS2=$(ls ${NAME}_R2*)
$bbmap/bbduk.sh in1=$i in2=$READS2 \
out1=${NAME}_R1_trimmed_bbduk.fastq.gz out2=${NAME}_R2_trimmed.fastq.gz \
ref=../ARTICv3_primers.fasta rcomp=f k=20 ktrim=r qtrim=rl trimq=3 maq=15 minlen=50 tbo
done

### TRIMMING UNPAIRED READS ###
for i in $(ls *ONT.fastq.gz);
do
NAME=${i/.fastq.gz/}
$bbmap/bbduk.sh in=$i out=${NAME}_trimmed.fastq.gz \
ref=../ARTICv3_primers.fasta rcomp=f k=20 ktrim=r qtrim=rl trimq=3 minlen=50 qin=33
done

for i in $(ls *IonTorrent.fastq.gz);
do
NAME=${i/.fastq.gz/}
$bbmap/bbduk.sh in=$i out=${NAME}_trimmed.fastq.gz \
ref=../ARTICv3_primers.fasta rcomp=f k=20 ktrim=r qtrim=rl trimq=3 minlen=50 qin=33
done



### ALIGNING LONG READS ###
minimap2 -d ref.mmi $REF # indexing reference
for i in $(ls *ONT_trimmed.fastq.gz);
do
NAME=${i/_trimmed.fastq.gz/}
minimap2 -ax map-ont $REF $i | samtools view -Sb - | samtools sort - > ${NAME}.bam 
done

### ALIGNING SINGLE-END SHORT READS (ION TORRENT) ###
for i in $(ls *IonTorrent_trimmed.fastq.gz);
do
NAME=${i/_trimmed.fastq.gz/}
minimap2 -ax sr $REF $i | samtools view -Sb - | samtools sort - > ${NAME}.bam 
done

### ALIGNING PAIRED-END SHORT READS ###
for i in $(ls *Illumina_R1_trimmed.fastq.gz);
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
# calling with bbmap callvariants.sh
ls -1 *.bam > bam-list
bamlist=./bam-list
DATE=$(date +'%Y%m%d')

$bbmap/callvariants.sh list=$bamlist out=alltimepoints_variants_${DATE}.vcf \
ref=$REF samstreamer=t ss=4 multisample=t clearfilters \
ploidy=1 mincov=0 overwrite=t

mv individual* ..
cd ..
bgzip -d individual*.vcf.gz
