#!/bin/bash

if [ $1 == align ] ; then
    mkdir -p ~/data/align
    for i in jsc bc bcbl ;
    do
	bowtie2-build -q ~/data/ref/$i.fasta ~/data/ref/$i
	bowtie2 -p 16 -x ~/data/ref/$i -1 ~/data/reads/trimmed/${i}_fwd_paired.fq.gz -2 ~/data/reads/trimmed/${i}_rev_paired.fq.gz | samtools view -bS - | samtools sort -o ~/data/align/$i.sorted.bam
	samtools index ~/data/align/$i.sorted.bam
    done
fi
