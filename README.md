# Visualizing Mutations in a Prolonged Infection with SARS-CoV-2

This repo contains the scripts used to produce the visualizations in the publication _Emergence of a globally unique SARS-CoV-2 Spike E484T mutation in a persistently infected immunocompromised individual._. These visualizations include:
1. __Figure 1A__ - Ct values through the course of a prolonged SARS-CoV-2 infection lasting 433+ days
2. __Figure 2A__ - Single-nucleotide variants (SNVs) in the consensus sequence of the patient's virus, shown at 8 time points. These SNVs are shown in their position in the SARS-CoV-2 Genome.
3. __Figure 2C__ - Results of a neutralization assay showing that on post-diagnosis day 297, the patient's virus was able to escape the Eli Lilly monoclonal antibody Bamlanivumab. This is the same date when the globally unique Spike E484T mutation was first detected.
4. __Supplemental Figure 1__ - Intrahost Single Nucleotide Variants (iSNVs) that reach consensus frequency (≥ 0.5 in sequencing reads) on post-diagnosis day 297, like Spike E484T did. Notably, these iSNVs all increased in frequency after the monoclonal therapy.
5. __Supplemental Figure 2__ - Count of genetic differences from Wuhan-1 in the patient's virus compared to genetic differences in 5,000 GISAID sequences, sampled evenly across the globe during the time course of the prolonged infection documented here.

All R scripts come with instructions on their usage, and are meant to be used within this repo as a working directory. To reproduce our results, we recommend downloading this repo as a ZIP, decompressing, and using it as a working directory.

Planned updates to this repo include:
* Switching to docker images of BBTools and Nextalign - DONE
* More instructions on how to reproduce each analysis/visualization, perhaps with a workflow manager
* Open access to data used in our study
