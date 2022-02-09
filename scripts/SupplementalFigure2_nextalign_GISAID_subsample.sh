#!/bin/bash


WD=${1:-$(locate prolonged_infection_workflow_wrapper_dockerized.sh | grep $(whoami) | sed 's/prolonged_infection_workflow_wrapper_dockerized.sh//')}
input="data/b12_enriched_global/b12_enriched_global_subsampled_sequences.fasta"
output="data/b12_enriched_global"
prefix="b12_enriched_global"

NEXTSTRAIN_DOCKER=(docker run -it --user $(id -u):$(id -g) -v $(pwd)/:/scratch -w /scratch nextstrain/base:build-20211210T215250Z)

cd $WD
$NEXTSTRAIN_DOCKER nextalign --sequences $input --reference data/reference.fasta --output-basename $prefix --include-reference -d $output
