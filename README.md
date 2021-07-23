# junco_flexibility_scripts

This repository includes scripts to replicate the analyses in:

**Stager M, Senner NR, Swanson DL, Carling MD, Eddy DK, Grieves TJ, and Cheviron ZA. 2021. Temperature heterogeneity correlates with intraspecific variation in physiological flexibility in a small endotherm. Nature Communications, accepted.**


Included are:

* field_data.R  
  * *R script for running analyses on phenotypic data from 335 individuals assayed in the field*
  * *This script draws on Sheet1 in the Source Data file*
  * *This script takes minutes to run on a 2012 MacBook Pro*

* RADSeq_processing 
  * *Commands for processing RAD-Seq data from 192 Junco sp. (demultiplexing w/ STACKS v.2.1, map w/ bwa, variant calling w/ STACKS, filtering w/ vcftools)*
  * *Associated bash commands in genome_map.sh, convertsam.sh, and sort_bam.sh*
  * *Raw reads used here are available from the NCBI Sequence Read Archive (PRJNA678344)*
  * *Barcodes provided in Junco_barcodes_SAMN16793657.txt and Junco_barcodes_SAMN16793658.txt (also in Sheet2 in the Source Data file)*
  * *Pop_map186.txt is a population map used by STACKS v.2.1 for 186 individuals (excluding 6 individuals that failed to sequence)*
  * *This script can take a couple days to run on a server with 32 cores, ~386 GB RAM*

* RDA.R 
  * *R script to run redundancy analysis on population genetic data*
  * *This script uses a plink file containing 21971 SNPs generated with above RADSeq_processing commands*
  * *This script draws on Sheet2 in the Source Data file*
  * *This script takes minutes to run on a 2012 MacBook Pro*
  
* pop_flex.R 
  * *R script for running analyses on phenotypic data from 95 acclimated individuals*
  * *This script draws on Sheet3 in the Source Data file*
  * *This script takes minutes to run on a 2012 MacBook Pro*

* RADSeq_PE_processing 
  * *Commands for processing paired-end RAD-Seq data from 95 individuals used in acclimation study*
  * *Associated bash commands in genome_map_pe.sh, convertsam.sh, sort_bam.sh, index_bam.sh, and merge_bam.sh*
  * *Raw reads used here are available from the NCBI Sequence Read Archive (PRJNA678344)*
  * *Barcodes provided in Junco_barcodes_SAMN16793659.txt (also in Sheet3 in the Source Data file)*
  * *Popmap_acc.txt is a population map used by STACKS v.2.1 for 89 individuals (excluding 6 individuals that failed to sequence); this file can be split into population-specific files for use in vcftools when calculating pairwise FST*
   * *This script can take a couple days to run on a server with 32 cores, ~386 GB RAM*
