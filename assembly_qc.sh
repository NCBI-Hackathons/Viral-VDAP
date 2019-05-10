#!/bin/bash

if [ $1 == contigs ] ; then
    echo num_contigs,n50,longest,shortest,total
    for i in jsc bc bcbl ;
    do
	python ~/Code/utils/qc/asm_assess.py -i ~/data/spades/$i/contigs.fasta
    done
fi

if [ $1 == scaffold ] ; then
    echo num_contigs,n50,longest,shortest,total
    for i in jsc bc bcbl ;
    do
	python ~/Code/utils/qc/asm_assess.py -i ~/data/medusa/$i/$i.contigs.fastaScaffold.fasta
    done
fi

if [ $1 == trim_mum ] ; then
    ##trim off ends of contigs using mummer info
    for i in jsc bc bcbl ;
    do
	nucmer -p ~/data/mummer/$i ~/data/ref/gk18.fasta ~/data/medusa/$i/$i.contigs.fastaScaffold.fasta 
	dnadiff -p ~/data/mummer/$i ~/data/ref/gk18.fasta ~/data/medusa/$i/$i.contigs.fastaScaffold.fasta
    done
fi


    
