### FIRST, run the installation on the miniconda3 installation here: https://docs.conda.io/en/latest/miniconda.html

conda update -n base conda
conda install -n base -c conda-forge mamba
mamba create -n nextstrain -c conda-forge -c bioconda \ # this creates an environment called nextstrain
  augur auspice nextstrain-cli nextalign snakemake awscli git pip # all this good stuff goes in the nextstrain environment

# to turn on the environment:
conda activate nextstrain
nextstrain check-setup --set-default # run this the first time you activate nextstrain

# finally, get all the coronavirus presets/defaults from github. Make sue you're aware of where the "ncov" folder goes:
git clone https://github.com/nextstrain/ncov.git
cd ncov
