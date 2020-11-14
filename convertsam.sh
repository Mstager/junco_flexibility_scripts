#!/bin/bash
SAMPLEID=`(ls ./sam/*.sam | sed 's/.sam//g')` #name sample ID as filename of raw reads without .fastq extension


INPUT_DIR=./sam #the input dir should have a /reads /trim and /sam folder



for i in $SAMPLEID
do echo "samtools for $i"
samtools view -S -b ${INPUT_DIR}/${i}.sam > ./bam/${i}.bam
done

exit 0;
