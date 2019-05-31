#!/bin/bash

if [ $1 == rmdup ] ; then
    ##remove duplicates
    for i in jsc bc bcbl ;
    do
	(samtools sort -n -o ~/data/align/$i.namesort.bam ~/data/align/$i.sorted.bam
	samtools fixmate -r -m -O bam ~/data/align/$i.namesort.bam ~/data/align/$i.rmd.bam 
	samtools sort -o ~/data/align/$i.rmd.sorted.bam ~/data/align/$i.rmd.bam
	samtools markdup -r -s -S ~/data/align/$i.rmd.sorted.bam ~/data/align/$i.nodup.bam
	samtools sort -o ~/data/align/$i.nodup.sorted.bam ~/data/align/$i.nodup.bam
	samtools index ~/data/align/$i.nodup.sorted.bam )&
    done
fi
