#!/bin/bash


if [ $1 == spades ] ; then
    mkdir -p ~/data/spades
    for i in jsc bc bcbl ;
    do
        spades.py -p 16 -x ~/data/reads/trimmed/${i}_fwd_paired.fq.gz -2 ~/data/reads/trimmed/${i}_rev_paired.fq.gz 

    done
fi

