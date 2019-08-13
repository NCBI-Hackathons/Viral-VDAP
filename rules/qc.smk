# path to reference genome
ref_fasta = config["ref"]["genome"]

##### Rule #####

rule fastqc:
    """Generate FastQC reports"""
    input:
        r1 = "data/reads/{sample}_r1.fq.gz",
        r2 = "data/reads/{sample}_r2.fq.gz" 
    output: 
        "data/qc_reports/{sample}/fastqc/{sample}_r1_fastqc.html",
        "data/qc_reports/{sample}/fastqc/{sample}_r2_fastqc.html"
    params: 
        "data/qc_reports/{sample}/fastqc"
    log: 
        "data/qc_reports/{sample}/fastqc/logs/{sample}.log"
    message: "Generating FastQC reports." 
    shell:
        "fastqc {input.r1} {input.r2} --outdir {params:} 2> {log}"

##### Rule #####

rule nucmer:
    """Generate delta file to be used to make plots with MUMmer"""
    input:
        ref = ref_fasta,
        scaffold = "data/assemblies/{sample}/scaffolds.fasta"
    output:
        "data/qc_reports/{sample}/mummerplot/{sample}.delta"
    params:
        prefix = "--prefix data/qc_reports/{sample}/mummerplot/{sample}"
    log:
        "data/qc_reports/{sample}/mummerplot/logs/nucmer.log"
    message: "Generating delta file with MUMmer"
    shell:
        "nucmer {params.prefix} {input.ref} {input.scaffold} 2> {log}"

##### Rule #####

rule mummerplot:
    """Generate plots with MUMmer to assess SPAdes assemblies"""
    input:
        ref = ref_fasta,
        scaffold = "data/assemblies/{sample}/scaffolds.fasta",
        delta = "data/qc_reports/{sample}/mummerplot/{sample}.delta"
    output:
        "data/qc_reports/{sample}/mummerplot/{sample}.png"        
    params: 
        prefix = "--prefix data/qc_reports/{sample}/mummerplot/{sample}", 
        plot = "--layout --filter --png"
    log:
        "data/qc_reports/{sample}/mummerplot/logs/mummerplot.log"
    message: "Generating plots aligning scaffolds to reference genome with MUMmer."
    shell:
        "mummerplot {params.prefix} {params.plot} -R {input.ref} -Q {input.scaffold} {input.delta} 2> {log}"  
