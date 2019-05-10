# dsVirus-variant-discovery-and-annotation-pipeline

## Team

Elena Maria Cornejo Castro - team lead  
Eneida Hatcher - writer  
Sara Jones  
Sasha Mushegian - writer  
Yunfan Fan  
Rashmi Naidu

# Introduction

Whole-genome sequencing of pathogenic viruses has the potential to improve surveillance, classification of disease subtypes, and association of viral subtypes with disease mechanisms. In order for viral genomic data to be universally interpretable and comparable, however, best practices must be established for quality control and variant calling. This is especially challenging for gamma-herpesvirus genomes, which are relatively large (>170 kb) and contain integrated human genes. Whole-genome analysis of these viruses is typically done on an ad-hoc basis by researchers working in isolation, making it difficult to know what types of comparative analyses are possible.

We have constructed a pipeline using freely available tools for quality control, alignment and SNP calling of double-stranded DNA virus paired-end short reads with the aim of providing researchers with interoperable consensus sequences and variant lists to be used for downstream analyses.    

## Is this the right pipeline for you?

**Your virus**
- Very low diversity virus populations, such as gammaherpesviruses. This tool is not suitable for viruses believed to exist as quasispecies (e.g., most RNA viruses), or for calling minority variants.
- Tested on data from enriched samples which have a higher proportion of viral reads compared to metagenomics or non-cultured samples.
- Input data should consist of paired-end reads from targeted sequencing. There is not a step to search for & remove host reads, since there should be very few host reads after targeted sequencing. Host reads are expected to be in a low enough concentration that they will be excluded when the consensus is built.

**What you’ll need:**
- FASTQ files of your NGS sequencing results
- FASTA files of reference sequence
-- Assuming you know the species you sampled,
-- Go to https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/find-data/virus and select “Search by Virus”
-- Begin typing the name of your virus. You can use taxonomic groups (e.g., Human gammaherpesvirus 4) or common names (e.g., Kaposi's sarcoma-associated herpesvirus, don’t worry, it’s an autofill, you don’t have to type the whole thing). 
-- On the filter panel on the left, click “Nucleotide Sequence Type,” then check “RefSeq.” Select the sequence you want, then download the FASTA file.

## Overview of pipeline steps
- Quality filtering, trimming, and minimum length filtering (Trimmomatic)
- De novo sequence assembly (SPAdes) -> Align scaffolds to reference and condense aligned scaffolds into a consensus/draft genome (Medusa)
- (Alternative to de novo assembly: Alignment to reference sequence (Bowtie2)-> Sequence deduplication (Samtools) -> calling variants and making consensus sequence (Samtools))
**Note:** On test data, de novo assembly produced a more complete assembly that better reproduced the corresponding published genome sequence than reference-based alignment. Steps for reference-based alignment are provided for the curious.

- Variant calling: 


## Software citations, versions and parameters
**Trimmomatic v.0.39** 
- Bolger, A. M., Lohse, M., & Usadel, B. (2014). Trimmomatic: A flexible trimmer for Illumina Sequence Data. Bioinformatics, btu170

Parameters: `ILLUMINACLIP:adapters.fa:2:30:10 LEADING:3 TRAILING:3 AVGQUAL:30 MINLEN:50`
- This line in the 1.trim.sh specifies that reads are to have a minimum average quality score of 30, low quality (<3) leading and trailing bases trimmed, and a minimum length after trimming of 50. All other settings as default. 

**SPAdes v3.13.0**

Parameters: `-k 21,33,55,77 -t 10 --only-assembler --careful`
- The line listed in sbatch.sh specifies that SPAdes should run in assembly module only (--only-assembler) and applies --careful to try to reduce the number of mismatches and short indels. The k parameter refers to the k-mer sizes. We used a range of sizes from 21 to 77. The t parameter refers to the number of threads to run the software. 

**Medusa v1.6**  
`All default settings`

For reference-based alignment (alternative):

**Bowtie2 v2.3.5.1**  
- Langmead B, Salzberg S. Fast gapped-read alignment with Bowtie 2. Nature Methods. 2012, 9:357-359.

`All default settings`  

**Samtools**  
- `3.consensus.sh` (`All default settings`)
- deduplication
- consensus generation: vcfutils and consensuscall-c to create a consesus sequence (FASTQ) from bowtie-produced alignments, then convert to FASTA


# How to use <this software>

# Installation options:

### Docker

The Docker image contains <this software> as well as a webserver and FTP server in case you want to deploy the FTP server. It does also contain a web server for testing the <this software> main website (but should only be used for debug purposes).

1. `docker pull ncbihackathons/<this software>` command to pull the image from the DockerHub
2. `docker run ncbihackathons/<this software>` Run the docker image from the master shell script
3. Edit the configuration files as below

# Testing

### DockerFile

<this software> comes with a Dockerfile which can be used to build the Docker image.

  1. `git clone https://github.com/NCBI-Hackathons/<this software>.git`
  2. `cd server`
  3. `docker build --rm -t <this software>/<this software> .`
  4. `docker run -t -i <this software>/<this software>`
