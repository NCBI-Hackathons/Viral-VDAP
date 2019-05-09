# dsVirus-variant-discovery-and-annotation-pipeline

## Team

Elena Maria Cornejo Castro - team lead  
Eneida Hatcher - writer  
Sara Jones  
Sasha Mushegian - writer  
Yunfan Fan  
Rashmi Naidu

## Introduction

Whole-genome sequencing of pathogenic viruses has the potential to improve surveillance, classification of disease subtypes, and association of viral subtypes with disease mechanisms. In order for viral genomic data to be universally interpretable and comparable, however, best practices must be established for quality control and variant calling. This is especially challenging for gamma-herpesvirus genomes, which are relatively large (>170 kb) and contain integrated human genes. Whole-genome analysis of these viruses is typically done on an ad-hoc basis by researchers working in isolation, making it difficult to know what types of comparative analyses are possible.

We have constructed a pipeline using freely available tools for quality control, alignment and SNP calling of double-stranded DNA virus paired-end short reads with the aim of providing researchers with interoperable consensus sequences and variant lists to be used for downstream analyses.    

# What is <this software>?

Overview Diagram

# How to use <this software>

# Software Workflow Diagram

# File structure diagram 
#### _Define paths, variable names, etc_

# Installation options:

We provide two options for installing <this software>: Docker or directly from Github.

### Docker

The Docker image contains <this software> as well as a webserver and FTP server in case you want to deploy the FTP server. It does also contain a web server for testing the <this software> main website (but should only be used for debug purposes).

1. `docker pull ncbihackathons/<this software>` command to pull the image from the DockerHub
2. `docker run ncbihackathons/<this software>` Run the docker image from the master shell script
3. Edit the configuration files as below

### Installing <this software> from Github

1. `git clone https://github.com/NCBI-Hackathons/<this software>.git`
2. Edit the configuration files as below
3. `sh server/<this software>.sh` to test
4. Add cron job as required (to execute <this software>.sh script)

### Configuration

```Examples here```

# Testing

We tested four different tools with <this software>. They can be found in [server/tools/](server/tools/) . 

# Additional Functionality

### DockerFile

<this software> comes with a Dockerfile which can be used to build the Docker image.

  1. `git clone https://github.com/NCBI-Hackathons/<this software>.git`
  2. `cd server`
  3. `docker build --rm -t <this software>/<this software> .`
  4. `docker run -t -i <this software>/<this software>`
  
### Website

There is also a Docker image for hosting the main website. This should only be used for debug purposes.

  1. `git clone https://github.com/NCBI-Hackathons/<this software>.git`
  2. `cd Website`
  3. `docker build --rm -t <this software>/website .`
  4. `docker run -t -i <this software>/website`
