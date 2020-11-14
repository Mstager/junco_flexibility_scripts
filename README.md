# junco_flexibility_scripts

This repository includes scripts to replicate the analyses in:

**Stager M, Senner NR, Swanson DL, Carling MD, Grieves TJ, and Cheviron ZA. 2020. Temperature heterogeneity correlates with intraspecific variation in physiological flexibility in a small endotherm. bioRxiv.**


Included are:

* junco_insitu.R  
  * *R script for running analyses on in situ data*
  * *This script draws on TableS1.csv from the Supplementary Materials*

* RADSeq_processing 
  * *Commands for processing RAD-Seq data from 192 Junco sp. (demultiplexing w/ STACKS, map w/ bwa, variant calling w/ STACKS, filtering w/ vcftools)*
  * *Associated bash commands in genome_map.sh, convertsam.sh, and sort_bam.sh
  * *Raw reads used here have been submitted to the NCBI Sequence Read Archive (PRJNA678344)*
  * *Barcodes provided in Junco_barcodes_lane1.txt and Junco_barcodes_lane2.txt (also in TableS2.csv in the Supp. Mat.)
  * *Pop_map186.txt is a population map used by STACKS for 186 individuals (6 individuals that failed to sequence not included)*

* RDA.R 
  * *R script to run redundancy analysis on population genetic data 
  * *This script uses Junco_SNPs.raw, a plink file containing 21971 SNPs generated with above RADSeq_processing commands
  * *This script draws on TableS2.csv from the Supplementary Materials*
  
* pop_flex.R 
  * *R script for running analyses on acclimated individuals*
  * *This script draws on TableS3.csv from the Supplementary Materials*

* RADSeq_PE_processing 
  * *Commands for processing RAD-Seq data from 95 individuals used in acclimation study*
  * *Associated bash commands in genome_map_pe.sh, convertsam.sh, and sort_bam.sh
  * *Raw reads used here have been submitted to the NCBI Sequence Read Archive (PRJNA678344)*
