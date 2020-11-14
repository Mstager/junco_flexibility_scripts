#!/bin/bash
SAMPLEID=`(ls *.rem.1.fq.gz | sed 's/.rem.1.fq.gz//g')` #name sample ID as filename of raw reads without .fastq extension


REF=/data/raw/maria/Junco_ref/GCA_003829775.1_ASM382977v1_genomic.fna #the location of the reference genome
INPUT_DIR=./ #the input dir should have a /reads /trim and /sam folder



for i in $SAMPLEID
do echo "bwa for $i"
bwa mem -t 20 ${REF} ${INPUT_DIR}/${i}.1.fq.gz ${INPUT_DIR}/${i}.2.fq.gz  > ${INPUT_DIR}/sam/${i}.paired.sam
done

exit 0;
