# Viral VDAP
Viral alignment, Variant Discovery and Annotation Pipeline

## Hackathon Team 
### First Women-led Biodata Science Hackathon <br> NIH campus, Bethesda, Maryland May 8-10

Elena Maria Cornejo Castro - team lead  
Eneida Hatcher - writer  
Sara Jones  
Sasha Mushegian - writer  
Yunfan Fan - sysadmin   
Rashmi Naidu

Allissa Dillman - hackathon organizer

# Introduction

Whole-genome sequencing of pathogenic viruses has the potential to improve surveillance, classification of disease subtypes, and association of viral subtypes with disease mechanisms. In order for viral genomic data to be universally interpretable and comparable, however, best practices must be established for quality control and variant calling. This is especially challenging for gamma-herpesvirus genomes, which are relatively large (>170 kb) and contain integrated human genes. Whole-genome analysis of these viruses is typically done on an ad-hoc basis by researchers working in isolation, making it difficult to know what types of comparative analyses are possible.

We have constructed a pipeline using freely available tools for quality control, alignment and SNP calling of double-stranded DNA virus paired-end short reads with the aim of providing researchers with interoperable consensus sequences and variant lists to be used for downstream analyses.    

# Is this the right pipeline for you?

**Your virus**
- Very low diversity virus populations, such as gammaherpesviruses. This tool is not suitable for viruses believed to exist as quasispecies (e.g., most RNA viruses), or for calling minority variants.
- Tested on data from enriched samples which have a higher proportion of viral reads compared to metagenomics or non-cultured samples.
- Input data should consist of paired-end reads from targeted sequencing. There is not a step to search for & remove host reads, since there should be very few host reads after targeted sequencing. Host reads are expected to be in a low enough concentration that they will be excluded when the consensus is built.

**What you’ll need:**
- FASTQ files of your NGS sequencing results
- FASTA files of reference sequence. If you don't have your own reference genome, read on. Assuming you know the species you sampled, 1) Go to https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/find-data/virus and select “Search by Virus" 2) Begin typing the name of your virus. You can use taxonomic groups (e.g., Human gammaherpesvirus...) or common names (e.g., Kaposi's sarcoma-associated herpesvirus... don’t worry, it’s an autofill, you don’t have to type the whole thing). 3) On the filter panel on the left, click “Nucleotide Sequence Type,” then check “RefSeq.” 4) Select the sequence you want, then download the FASTA file.

# Overview of pipeline steps
- QC: Quality filtering, trimming, and minimum length filtering (Trimmomatic)
- Assembly: De novo sequence assembly (SPAdes) -> Align scaffolds to reference and condense aligned scaffolds into a consensus/draft genome (Medusa)
- Variant discovery: Compare assembly to reference genome to extract SNPs (Parsnp) _[testing in progress]_

- **Note:** On test data, de novo assembly produced a more complete assembly that better reproduced the corresponding published genome sequence than reference-based alignment. For the curious, reference-based alignment can be carried out as follows (more details below): Alignment to reference sequence (Bowtie2)-> Sequence deduplication (Samtools) -> calling variants and making consensus sequence (Samtools))

## Software citations, versions and parameters
**Trimmomatic v.0.39** 
- Bolger, A. M., Lohse, M., & Usadel, B. (2014). Trimmomatic: A flexible trimmer for Illumina Sequence Data. Bioinformatics, btu170

Parameters: `ILLUMINACLIP:adapters.fa:2:30:10 LEADING:3 TRAILING:3 AVGQUAL:30 MINLEN:50`
- This line in `1.trim.sh` specifies that reads are to have a minimum average quality score of 30, low quality (<3) leading and trailing bases trimmed, and a minimum length after trimming of 50. All other settings as default. 

**SPAdes v3.13.0**
-  Bankevich A., Nurk S., Antipov D., Gurevich A., Dvorkin M., Kulikov A. S., Lesin V., Nikolenko S., Pham S., Prjibelski A., Pyshkin A., Sirotkin A., Vyahhi N., Tesler G., Alekseyev M. A., Pevzner P. A. SPAdes: A New Genome Assembly Algorithm and Its Applications to Single-Cell Sequencing.    Journal of Computational Biology, 2012

Parameters: `-k 21,33,55,77 -t 10 --only-assembler --careful`
- The line listed in `sbatch.sh` specifies that SPAdes should run in assembly module only (--only-assembler) and applies --careful to try to reduce the number of mismatches and short indels. The k parameter refers to the k-mer sizes. We used a range of sizes from 21 to 77. The t parameter refers to the number of threads to run the software. 

**Medusa v1.6**  
`scaffold.sh` (`All default settings`)
Downstream steps use the largest scaffold produced from alignment.

**Parsnp v1.2** 
- Treangen TJ*, Ondov BD*, Koren S, Phillippy AM:
     Rapid Core-Genome Alignment and Visualization for Thousands of Microbial Genomes.
     bioRxiv (2014). doi: http://dx.doi.org/10.1101/007351

Parameters: `All default settings`

**Snakemake v4.3.1**
- Köster, Johannes and Rahmann, Sven. “Snakemake - A scalable bioinformatics workflow engine”. Bioinformatics 2012.

Parameters: `All default settings`

To check for errors in Snakefile or config.yaml:

`snakemake -np --configfile config.yaml`

To run snakemake:

`snakemake --configfile config.yaml`

For reference-based alignment (alternative):

_Bowtie2 v2.3.5.1_  
- Langmead B, Salzberg S. Fast gapped-read alignment with Bowtie 2. Nature Methods. 2012, 9:357-359.

`All default settings`  

_Samtools_  
- `3.consensus.sh` (`All default settings`)
- deduplication
- consensus generation: vcfutils and consensuscall-c to create a consesus sequence (FASTQ) from bowtie-produced alignments.

_SeqTK_
- (`All default settings`)
- converting FASTQ consensus sequence to FASTA format


# How to use <this software>

# Installation options:

### Docker

The Docker image contains <this software> as well as a webserver and FTP server in case you want to deploy the FTP server. It does also contain a web server for testing the <this software> main website (but should only be used for debug purposes).

1. `docker pull ncbihackathons/<this software>` command to pull the image from the DockerHub
2. `docker run ncbihackathons/<this software>` Run the docker image from the master shell script
3. Edit the configuration files as below

### DockerFile

<this software> comes with a Dockerfile which can be used to build the Docker image.

  1. `git clone https://github.com/NCBI-Hackathons/<this software>.git`
  2. `cd server`
  3. `docker build --rm -t <this software>/<this software> .`
  4. `docker run -t -i <this software>/<this software>`

# Testing

We tested sequences from three different KSHV cell lines with our pipeline. These data were derived by the Wellcome Sanger Institute and are available in the European Nucleotide Archive [JSC-1 (SAMEA1709534), BC-1 (SAMEA1709542), BCBL-1 (SAMEA1709549)]. The reference genome used was the KSHV GK18 strain complete genome sequence (Accession number: AF148805). The genome for GK18 is currently the most extensively annotated KSHV sequence available, including gene and coding sequence, repeat regions, mRNA and PolyA features (Rezaee et al., 2006).

![Progressive Mauve alignment of JSC-1 de novo sequence.](https://github.com/NCBI-Hackathons/dsVirus-variant-discovery-and-annotation-pipeline/blob/master/images/JSC-1_Mauve.png)
Progressive Mauve alignment of JSC-1 de novo sequence.

The newly derived JSC-1 sequence was aligned to GK18 (AF148805) and the JSC-1 genome deposited in RefSeq (GQ994937). GK18 annotations are represent as blocks, with gene regions in white and known repeat regions in red.
As expected, the de novo assembly does not resolve the repeat regions.



