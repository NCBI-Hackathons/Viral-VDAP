### Set of rules to create consensus FASTA file that will be used for calling variants against reference genome ###

rule samtools_index:
    """Index BAM file with SAMtools"""
    input: 
        "data/dedup_reads/{sample}_dedup.bam"
    output:
        "data/dedup_reads/{sample}_dedup.bam.bai"
    message: "Indexing BAM file with SAMtools."
    shell:
        "samtools index {input}"


rule create_consensus:
    """Call variants on BAM file to create consensus FASTA with SAMtools, BCFtools, and Seqtk"""    
    input:
        fa = "data/final_scaffolds/{sample}/{sample}_final_scaffold.fasta",
        bam = "data/dedup_reads/{sample}_dedup.bam",
        bai = "data/dedup_reads/{sample}_dedup.bam.bai"
    output:
        "data/consensus_fasta/{sample}/{sample}_consensus.fasta"
    message: "Generating consensus fasta with SAMtools, BCFtools, Seqtk."""
    shell: 
        "bcftools mpileup -f {input.fa} {input.bam} | bcftools call -cM - | vcfutils.pl vcf2fq | seqtk seq -A - > {output}" 


rule variant_calling:
    """Call variants on consensus FASTA file with Parsnp"""
    input: 
        ref = ref,
        sample = "data/consensus_fasta/{sample}/{sample}_consensus.fasta"
    output: 
        "data/vcf/{sample}/parsnp.ggr",
        "data/vcf/{sample}/parsnp.tree",
        "data/vcf/{sample}/parsnp.xmfa",
        "data/vcf/{sample}/parsnpAligner.log",
        "data/vcf/{sample}/parsnpAligner.ini"  
    threads: num_threads 
    log:
        "data/vcf/{sample}/parsnp.log"
    message: "Executing variant calling with Parsnp."
    params: 
        sample_dir = "data/consensus_fasta/{sample}/",
        output_dir = "data/vcf/{sample}/"
    shell: 
        "./Parsnp-Linux64-v1.2/parsnp -r {input.ref} -d {params.sample_dir} -p {threads} -c -o {params.output_dir} 2> {log}"

rule generate_vcf:
    """Generate VCF from parsnp.ggr file with HarvestTools"""
    input:
        ggr ="data/vcf/{sample}/parsnp.ggr",
        xmfa = "data/vcf/{sample}/parsnp.xmfa",
        tree = "data/vcf/{sample}/parsnp.tree",
        ini = "data/vcf/{sample}/parsnpAligner.ini",
        log = "data/vcf/{sample}/parsnpAligner.log"
    output:
        "data/vcf/{sample}/{sample}.vcf"
    message: "Executing HarvestTools to generate VCF file."
    shell:
        "harvesttools-Linux64-v1.2/harvesttools -i {input.ggr} -V {output}"

