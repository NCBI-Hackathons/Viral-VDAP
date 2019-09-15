rule denovo_assembly:
    """Denovo assembly with SPAdes."""
    input:
        r1 = "data/trimmed_reads/{sample}_forward_paired.fq.gz",
        r2 = "data/trimmed_reads/{sample}_reverse_paired.fq.gz"
    output: 
        "data/assemblies/{sample}/scaffolds.fasta"
    params: 
        prefix = "data/assemblies/{sample}",
        extra = add_ons.strip() # remove leading and trailing white space
    threads: num_threads
    message: "Executing denovo assembly using SPAdes."
    # checks if params is empty string    
    run:
        if not params.extra: 
            shell("SPAdes-3.5.0-Linux/bin/spades.py -t {threads} --careful -1 {input.r1} -2 {input.r2} -o {params.prefix}")
        else:
            shell("SPAdes-3.5.0-Linux/bin/spades.py {params.extra} -t {threads} --careful -1 {input.r1} -2 {input.r2} -o {params.prefix}") 
