rule trim_and_filter_reads:
    """
    Trims and filters reads from paired-end FASTQ files using Trimmomatic.
    """
    input:
        unpack(get_fastq)
    output:
        temp("trimmed/{sample}_forward_paired.fq.gz"),
        temp("trimmed/{sample}_forward_unpaired.fq.gz"),
        temp("trimmed/{sample}_reverse_paired.fq.gz"),
        temp("trimmed/{sample}_reverse_unpaired.fq.gz")
    params:
        config["params"]["trimmomatic"]["options"]
    log:
        "logs/trimmomatic/{sample}-{unit}.log"
    threads: config["params"]["trimmomatic"]["threads"] 
    shell:
        "java -jar /path/of/Trimmomatic-0.39/trimmomatic-0.39.jar PE -t {threads} {input} {output} 2> {log}"
