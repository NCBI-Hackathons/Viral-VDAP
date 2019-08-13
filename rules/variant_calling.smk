rule samtools_index:
    """Index BAM file with SAMtools"""
    input: 
        "data/dedup_reads/{sample}_dedup.bam"
    output:
        "data/dedup_reads/{sample}_dedup.bam.bai"
    params:
        file = "{sample}_dedup.bam.bai",
        dir = "data/dedup_reads"
    message: "Indexing BAM file with SAMtools."
    shell:
        "samtools index {input} && mv {params.file} {params.dir}"


rule bcftools_call:
    """Call SNPS and INDELS to consensus fasta with BCFtools"""    
    input:
        fa = "data/consensus_fasta/{sample}/{sample}_consensus.fasta",
        bam = "data/dedup_reads/{sample}_dedup.bam",
        bai = "data/dedup_reads/{sample}_dedup.bam.bai"
    output:
        "data/vcf/{sample}_raw.vcf"
    message: "Calling variants and indels against consensus fasta with BCFtools."""
    shell: 
        "samtools mpileup -g -f {input.fa} -a -s {input.bam} | bcftools call -mM - > {output}" 
   
