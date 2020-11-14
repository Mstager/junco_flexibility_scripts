#!/bin/bash
SAMPLEID=`(ls *.rem.1.fq.gz | sed 's/.rem.1.fq.gz//g')` #name sample ID as filename 


INPUT_DIR=./bam_sorted #the input dir should have a /reads /trim and /sam folder



for i in $SAMPLEID
do echo "samtools for $i"
samtools merge ./merged_mapped_reads/${i}.bam ${INPUT_DIR}/${i}.paired.bam ${INPUT_DIR}/${i}.rem.1.bam ${INPUT_DIR}/${i}.rem.2.bam 
done

exit 0;


