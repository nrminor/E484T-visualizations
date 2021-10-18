# pip install snipit

WD="/Users/nicholasminor/Documents/informatics/E484T_paper"

cd $WD
software/nextalign --sequences data/alltimepoints.fasta --reference data/ref/reference.fasta --genemap data/genemap.gff --genes S --output-basename nextalign --include-reference -d data/nextalign_output/

snipit data/nextalign_output/nextalign.aligned.fasta -d visuals/ -l data/labels.csv --l-header 'name,label'
snipit data/nextalign_output/nextalign.gene.S.fasta -d visuals -o spike_mutations -l data/labels.csv --l-header 'name,label'
