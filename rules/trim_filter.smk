rule merge_forward_reads:
    """
    Merge forward reads of the same sample from different runs if more than one.
    Rename FASTQ files with sample name.
    Compress files if not already compressed. 
    """
    input:
        unpack(get_forward_fastqs)
    output:
        "reads/{sample}_r1.fq.gz"
    run:
        if len(input) > 1:
            shell("cat {input} | gzip > {output}")
        else:
            shell("mv {input} {output}")
               
rule merge_reverse_reads:
    """
    Merge reverse reads of the same sample from different runs if more than one.
    Rename FASTQ files with sample name.
    Compress files if not already compressed. 
    """
    input:
        unpack(get_reverse_fastqs)
    output:
        "reads/{sample}_r2.fq.gz"
    run:
        if len(input) > 1:
            shell("cat {input} | gzip > {output}")
        else:
            shell("mv {input} {output}")

rule trim_and_filter_reads:
    """
    Trims and filters reads from paired-end FASTQ files using Trimmomatic.
    """
    input:
        "reads/{sample}_r1.fq.gz",
        "reads/{sample}_r2.fq.gz"
    output:
        temp("trimmed/{sample}_forward_paired.fq.gz"),
        temp("trimmed/{sample}_forward_unpaired.fq.gz"),
        temp("trimmed/{sample}_reverse_paired.fq.gz"),
        temp("trimmed/{sample}_reverse_unpaired.fq.gz")
    params:
        trim = " ".join(config["params"]["trimmomatic"]["options"]),
        clip = config["params"]["trimmomatic"]["clip"].format([config]["adapters"])
    log:
        "logs/trimmomatic/{sample}.log"
    threads: config["params"]["trimmomatic"]["threads"] 
    shell:
        "trimmomatic PE -t {threads} {input} {output} {params.trim} {params.clip} 2> {log}"
