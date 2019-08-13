### Set of rules to map reads, sort, and remove duplicates ###

rule bwa_index:
    """Index consensus FASTA file with BWA"""
    input: 
       "data/consensus_fasta/{sample}/{sample}_consensus.fasta" 
    output:
       "data/consensus_fasta/{sample}/{sample}.amb",
       "data/consensus_fasta/{sample}/{sample}.ann",
       "data/consensus_fasta/{sample}/{sample}.bwt",
       "data/consensus_fasta/{sample}/{sample}.pac",
       "data/consensus_fasta/{sample}/{sample}.sa" 
    log: 
       "data/consensus_fasta/{sample}.log"
    params:
       prefix = "{sample}",
       directory = "data/consensus_fasta/{sample}"
    message: "Indexing consensus FASTA file with BWA."
    shell:
        "./bwa-0.7.17/bwa index -p {params.prefix} {input} && mv {params.prefix}* {params.directory}"

rule bwa_map:
    """Map reads to consensus FASTA, outputs BAM file"""
    input:
        r1 = "data/trimmed_reads/{sample}_forward_paired.fq.gz",
        r2 = "data/trimmed_reads/{sample}_reverse_paired.fq.gz",
        ref_amb = "data/consensus_fasta/{sample}/{sample}.amb",
        ref_ann = "data/consensus_fasta/{sample}/{sample}.ann",
        ref_bwt = "data/consensus_fasta/{sample}/{sample}.bwt", 
        ref_pac = "data/consensus_fasta/{sample}/{sample}.pac", 
        ref_sa = "data/consensus_fasta/{sample}/{sample}.sa"
    output:
        "data/mapped_reads/{sample}.bam"
    params: 
        ref_base = "data/consensus_fasta/{sample}/{sample}"
    threads: num_threads
    message: "Mapping reads to consensus FASTA file with BWA."
    log:
        "data/mapped_reads/logs/{sample}.log"
    shell: 
        "(./bwa-0.7.17/bwa mem -t {threads} {params.ref_base} {input.r1} {input.r2} | samtools view -Sb - > {output}) 2> {log}"

rule samtools_sort:
    """Sort BAM file with SAMtools"""
    input: 
        "data/mapped_reads/{sample}.bam"
    output: 
        "data/sorted_reads/{sample}_sorted.bam"
    message: "Sorting reads in BAM file with SAMtools."
    shell: 
        "samtools sort -O bam {input} > {output}"

rule remove_duplicates:
    """Remove duplicates"""
    input:
        "data/sorted_reads/{sample}_sorted.bam"
    output: 
        "data/dedup_reads/{sample}_dedup.bam"
    message: "Removing duplicates from BAM files with Picard."
    log:
        metrics = "data/dedup_reads/logs/{sample}_metrics.log",
        stdout = "data/dedup_reads/logs/{sample}.log"
    shell:
        "java -jar picard.jar MarkDuplicates I={input} O={output} REMOVE_DUPLICATES=true METRICS_FILE={log.metrics} 2> {log.stdout}"
