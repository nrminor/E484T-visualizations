#!/bin/bash

DATE=$(date +'%Y%m%d')

while getopts :d:s:m:p: flag
do
	case "${flag}" in
		d) outdir=${OPTARG};;
		:) outdir=data/b12_enriched_global ;;
		s) sequence=${OPTARG};;
		m) metadata=${OPTARG};;
		i) include=${OPTARG};;
		p) prefix=${OPTARG};;
		:) prefix=gisaid_${DATE} ;;
	esac
done

SCRIPTS=scripts/ncov_full/scripts


# reformatting gisaid fasta deflines 
docker run -it --user $(id -u):$(id -g) -v $(pwd)/:/scratch -w /scratch nextstrain/base:build-20211210T215250Z \
python3 $SCRIPTS/sanitize_sequences.py \
	--sequences $sequence \
	--strip-prefixes "hCoV-19/" \
	--output ${outdir}/${prefix}.fasta.gz && \
echo "FASTA deflines successfully reformatted"

# indexing gisaid fasta
docker run -it --user $(id -u):$(id -g) -v $(pwd)/:/scratch -w /scratch nextstrain/base:build-20211210T215250Z \
augur index \
	--sequences ${outdir}/${prefix}.fasta.gz \
	--output ${outdir}/${prefix}_index.tsv.gz && \
echo "FASTA successfully indexed"

# renaming column headers and strain names in gisaid metadata
docker run -it --user $(id -u):$(id -g) -v $(pwd)/:/scratch -w /scratch nextstrain/base:build-20211210T215250Z \
python3 $SCRIPTS/sanitize_metadata.py \
	--metadata ${metadata} \
	--parse-location-field Location \
	--rename-fields 'Virus name=strain' 'Accession ID=gisaid_epi_isl' 'Collection date=date' \
	--strip-prefixes "hCoV-19/" \
	--output ${outdir}/${prefix}_metadata.tsv.gz && \
echo "Metadata successfully reformatted"

# creating include list from 500 closely related B.1.2 samples from UShER
# cut -f 1 $(locate usher_500samples_metadata.tsv) > ${outdir}/include_B12.txt 
cut -f 1 $include > ${outdir}/include_B12.txt && \
echo "UShER output successfully found and stored in" ${outdir}
echo "The strain list is called include_B12.txt"

# filtering down strain list to our desired subsample
docker run -it --user $(id -u):$(id -g) -v $(pwd)/:/scratch -w /scratch nextstrain/base:build-20211210T215250Z \
augur filter \
	--metadata ${outdir}/${prefix}_metadata.tsv.gz \
	--min-date 2020-09-01 \
	--max-date 2022-02-01 \
	--exclude-ambiguous-dates-by any \
	--subsample-max-sequences 5000 \
	--group-by region year month \
	--include ${outdir}/include_B12.txt \
	--output-strains ${outdir}/strains_global.txt && \
echo "Strain list successfully filtered down and stored in" ${outdir}
echo "The expanded strain list is called strains_global.txt"


# filtering fasta and metadata to match the subsample we defined above
docker run -it --user $(id -u):$(id -g) -v $(pwd)/:/scratch -w /scratch nextstrain/base:build-20211210T215250Z \
augur filter \
	--metadata ${outdir}/${prefix}_metadata.tsv.gz \
	--sequence-index ${outdir}/${prefix}_index.tsv.gz \
	--sequences ${outdir}/${prefix}.fasta.gz \
	--exclude-all \
	--include ${outdir}/strains_global.txt \
	--output-metadata ${outdir}/${prefix}_subsampled_metadata.tsv.gz \
	--output-sequences ${outdir}/${prefix}_subsampled_sequences.fasta.gz && \
echo "The FASTA and metadata have now been successfully subsampled."
ls ${outdir}/*subsampled*
echo "These files are ready to be input into the nextstrain workflow as follows:"
echo "inputs:"
echo "  - name: subsampled-gisaid"
echo "	metadata: " ${outdir}/${prefix}_subsampled_metadata.tsv.gz
echo "	sequences: " ${outdir}/${prefix}_subsampled_sequences.fasta.gz

gunzip ${outdir}/${prefix}_subsampled_sequences.fasta.gz
gunzip ${outdir}/${prefix}_subsampled_metadata.tsv.gz

rm ${outdir}/${prefix}.fasta.gz
rm ${outdir}/${prefix}_metadata.tsv.gz
rm ${outdir}/${prefix}_index.tsv.gz


# mapping subsample fasta to reference sequence
mkdir ${outdir}/tmp
Rscript $SCRIPTS/fasta_sep.R $(pwd)/${outdir}
docker run -i --user $(id -u):$(id -g) -v $(pwd)/:/scratch -w /scratch quay.io/biocontainers/minimap2:2.24--h5bf99c6_0 \
/bin/bash scripts/minimap_gisaid_fastas.sh

# calling VCFs from GISAID subsample
docker run -i --user $(id -u):$(id -g) -v $(pwd)/:/scratch -w /scratch quay.io/biocontainers/bbmap:38.93--he522d1c_0 \
callvariants.sh \
list=${outdir}/tmp/sam_list.txt out=${outdir}/${prefix}_subsampled_sequences.vcf.gz \
ref=ref/reference.fasta samstreamer=t ss=4 multisample=t clearfilters \
ploidy=1 mincov=0 callsub=t calldel=f callins=f overwrite=t
rm -r ${outdir}/tmp
gunzip ${outdir}/*.vcf.gz
