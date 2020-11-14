#!/bin/bash
SAMPLEID=`(ls *.fq.gz | sed 's/.fq.gz//g')` #name sample ID as filename of raw reads without .fastq extension


REF=/data/raw/maria/Junco_ref/GCA_003829775.1_ASM382977v1_genomic.fna #the location of the reference genome
INPUT_DIR=./ #the input dir should have a /reads /trim and /sam folder



for i in $SAMPLEID
do echo "bwa for $i"
bwa mem -t 10 ${REF} ${INPUT_DIR}/${i}.fq.gz > ${INPUT_DIR}/sam/${i}.sam
done

exit 0;
