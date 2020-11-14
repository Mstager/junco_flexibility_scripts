#!/bin/bash
SAMPLEID=`(ls ./bam/*.bam | sed 's/.bam//g')` #name sample ID as filename of raw reads without .fastq extension


INPUT_DIR=./bam_sorted #the input dir should have a /reads /trim and /sam folder



for i in $SAMPLEID
do echo "samtools for $i"
samtools index ${INPUT_DIR}/${i}.bam  ${INPUT_DIR}/${i}.bai
done

exit 0;
