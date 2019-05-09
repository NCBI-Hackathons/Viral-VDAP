#!/bin/bash

if [ $1 == scaffold ] ; then
    mkdir -p ~/data/medusa
    for i in jsc bc bcbl ;
    do
	mkdir -p ~/data/medusa/$i
	cp -r ~/software/medusa/* ~/data/medusa/$i/
	cp ~/data/spades/$i/contigs.fasta ~/data/medusa/$i/$i.contigs.fasta
	mkdir -p ~/data/medusa/$i/ref
	cp ~/data/ref/gk18.fasta ~/data/medusa/$i/ref
	cd ~/data/medusa/$i
	java -jar ./medusa.jar -f ~/data/medusa/$i/ref -i ~/data/medusa/$i/$i.contigs.fasta -v
    done
fi    
	     
