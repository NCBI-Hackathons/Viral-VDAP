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
