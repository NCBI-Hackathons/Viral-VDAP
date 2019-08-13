import pandas as pd

# Config file:
configfile: "config.yaml"

# Read in samples.tsv file
samples = pd.read_table(config["samples"]).set_index("sample",drop=False)
               
# Path to adapters file
adapters = config["adapters"]

# Number of threads to use
num_threads = int(config["threads"])

##### Helper functions #####

def get_fastq(wildcards):
    """Get fastq files of given sample."""
    return samples.loc[wildcards.sample, ["fq1", "fq2"]].dropna()

##### Rule #####

rule trim_and_filter_reads:
    """Trims and filters reads from paired-end FASTQ files using Trimmomatic."""
    input:
        get_fastq
    output:
        "data/trimmed_reads/{sample}_forward_paired.fq.gz",
        "data/trimmed_reads/{sample}_forward_unpaired.fq.gz",
        "data/trimmed_reads/{sample}_reverse_paired.fq.gz",
        "data/trimmed_reads/{sample}_reverse_unpaired.fq.gz"
    params:
        trim = "LEADING:3 TRAILING:3 AVGQUAL:30 MINLEN:50",
        clip = "ILLUMINACLIP:{}:2:30:10".format(adapters)
    log:
        "data/trimmed_reads/{sample}.log"
    threads: num_threads
    message: "Executing trimming of reads using Trimmomatic." 
    shell:
        "java -jar Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads {threads} -phred33 -trimlog {log} {input} {output} {params.trim} {params.clip}"
