#!/bin/bash
SAMPLEID=`(ls ./bam/*.bam | sed 's/.bam//g')` #name sample ID as filename of raw reads without .fastq extension


INPUT_DIR=./bam #the input dir should have a /reads /trim and /sam folder



for i in $SAMPLEID
do echo "samtools for $i"
samtools sort ${INPUT_DIR}/${i}.bam > ./bam_sorted/${i}.bam
done

exit 0;
