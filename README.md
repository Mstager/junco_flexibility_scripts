# junco_flexibility_scripts

This repository includes scripts to replicate the analyses in:

**Stager M, Senner NR, Swanson DL, Carling MD, Grieves TJ, and Cheviron ZA. 2020. Temperature heterogeneity correlates with intraspecific variation in physiological flexibility in a small endotherm. bioRxiv.**


Included are:

* junco_insitu.R  
  * *R script for running analyses on in situ data*
  * *This script draws on TableS1.csv from the Supplementary Materials*

* RADSeq_processing 
  * *Commands for processing RAD-Seq data (demultiplexing w/ STACKS, map w/ bwa, variant calling w/ STACKS, filtering w/ vcftools)*
  * *Associated bash commands in genome_map.sh, convertsam.sh, and sort_bam.sh

* Pop_map186.txt 
  * *Population map used by STACKS for 186 individuals (6 individuals that failed to sequence not included)

* Junco_SNPs.raw 
  * *plink file containing 21971 SNPs generated with above RADSeq_processing commands

* RDA.R 
  * *R script to run redundancy analysis on population genetic data 
  * *This script uses Junco_SNPs.raw*
  * *This script draws on TableS2.csv from the Supplementary Materials*
  
* pop_flex.R 
  * *R script for running analyses on acclimated individuals*
  * *This script draws on TableS3.csv from the Supplementary Materials*





Associated data are available in the supplementary materials:

1. TableS1.csv # *A file containing data for individuals sampled in situ*

2. TableS2.csv # *A file containing data for population genetic samples from across the Junco genus*

3. TableS3.csv # *A file containing data for acclimated individuals*

Raw reads will be available from the NCBI Sequence Read Archive (PRJNA678344).
