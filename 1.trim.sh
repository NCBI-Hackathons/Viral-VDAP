#!/bin/bash

if [ $1 == trimmo ] ; then
    datadir=~/data/reads
    readsdir=$datadir/raw
    mkdir -p $datadir/trimmed
    for i in jsc bcbl bc ;
    do
	java -jar ~/software/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 16 -phred33 \
	     $readsdir/${i}_r1.fastq.gz $readsdir/${i}_r2.fastq.gz \
	     $datadir/trimmed/${i}_fwd_paired.fq.gz \
	     $datadir/trimmed/${i}_fwd_unpaired.fq.gz \
	     $datadir/trimmed/${i}_rev_paired.fq.gz \
	     $datadir/trimmed/${i}_rev_unpaired.fq.gz \
	     ILLUMINACLIP:adapters.fa:2:30:10 LEADING:3 TRAILING:3 AVGQUAL:30 MINLEN:50
    done
fi
    
