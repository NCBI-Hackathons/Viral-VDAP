#!/usr/bin/env python3 

configfile: "config.yaml"

"""
rule all:
    input:
"""

rule trim_filter:
    """
    Uses Trimmomatic to trim out adapters and leading and trailing bases of low quality (<3) from
    forward and reverse fastq files.
    Reads are removed if they have average quality score < 30 and length < 50.

    """
    input:
        R1 = config["R1"],
        R2 = config["R2"],
        adapters = config["adapters"]
    output:
        "trimmed/jsc_fwd_paired.fq.gz",
        "trimmed/jsc_fwd_unpaired.fq.gz",
        "trimmed/jsc_rev_paired.fq.gz",
        "trimmed/jsc_rev_unpaired.fq.gz"
    shell:
        "java -jar ~/software/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 8 -phred33 {input.R1} {input.R2} {output} ILLUMINACLIP:{input.adapters}:2:30:10 LEADING:3 TRAILING:3 AVGQUAL:30 MINLEN:50"

'''
rule denovo_assemble:
    """
    Uses Spades for denovo assembly of the forward and reverse reads (fastq).

    """
    input:
    output:
    threads: 8 
    shell:
    spades.py -k 21,33,55,77 -t 10 \
         --only-assembler --careful  \
    -1 ~/data/reads/trimmed/${i}_fwd_paired.fq.gz \
    -2 ~/data/reads/trimmed/${i}_rev_paired.fq.gz \
    -o ~/data/spades/$i


rule medusa_align:
    """
    Uses Medusa to align scaffolds (fasta) to the reference genome (fasta).
    """
    input:
    output:
    shell:
        java -jar ./medusa.jar -f ~/data/medusa/$i/ref -i ~/data/medusa/$i/$i.contigs.fasta -v
 

rule find_snps:
    """
    Uses parsnps to find snps between aligned scaffolds (fasta) and the reference genome (fasta).
    """
    input:
    output:
    shell:
        parsnp -p 16 -c -d ~/data/scaffolds -r ~/data/scaffolds/gk18.fasta -o ~/data/parsnp
harvesttools -i ~/data/parsnp/parsnp.ggr -V ~/data/parsnp/parsnp.vcf

'''
