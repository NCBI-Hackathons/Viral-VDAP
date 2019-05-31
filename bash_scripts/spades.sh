#!/bin/bash


mkdir -p ~/data/spades
for i in jsc bc bcbl ;
do
    mkdir -p ~/data/spades/$i
    spades.py -k 21,33,55,77 -t 10 \
	      --only-assembler --careful  \
    -1 ~/data/reads/trimmed/${i}_fwd_paired.fq.gz \
    -2 ~/data/reads/trimmed/${i}_rev_paired.fq.gz \
       -o ~/data/spades/$i 
done


