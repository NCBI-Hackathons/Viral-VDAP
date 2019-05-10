#!/bin/bash

parsnp -p 16 -c -d ~/data/scaffolds -r ~/data/scaffolds/gk18.fasta -o ~/data/parsnp 
harvesttools -i ~/data/parsnp/parsnp.ggr -V ~/data/parsnp/parsnp.vcf
