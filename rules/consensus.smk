rule extract_scaffold:
    """Extract longest sequence from scaffolds.fasta to generate FASTA file to be used for mapping reads against"""
    input:
        "data/scaffolds/{sample}_scaffolds.fasta"
    output:
        "data/final_scaffolds/{sample}/{sample}_final_scaffold.fasta"
    message: "Creating consensus FASTA file."
    run:
        id = get_id(input)
        shell("samtools faidx {input} {id} > {output}")
