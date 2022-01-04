WD="/Users/nicholasminor/Documents/informatics/E484T_paper" ; input="/Volumes/working_ssd/e484t_manuscript/data/b12_enriched_global/b12_enriched_global_subsampled_sequences.fasta" ; output="/Volumes/working_ssd/e484t_manuscript/data/b12_enriched_global" ; prefix="b12_enriched_global"

cd $WD
software/nextalign --sequences $input --reference data/ref/reference.fasta --output-basename $prefix --include-reference -d $output
