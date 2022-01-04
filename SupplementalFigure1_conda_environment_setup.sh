# Recipe for Conda environment used for variant-calling

conda update -n base conda
conda install -n base -c conda-forge miniconda3

conda create --name variant-calling python=3.5 

conda activate variant-calling
conda install -c bioconda vcftools
conda install -c bioconda bcftools
conda install -c bioconda vcflib
conda install -c bioconda samtools
conda install -c bioconda bwa
conda install -c bioconda minimap2
conda install -c bioconda hisat2
conda install -c bioconda trimmomatic
conda install -c bioconda bedtools
