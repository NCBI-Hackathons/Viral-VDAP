rule refine_assembly:
    """Refine assembly using Medusa"""
    input: 
        ref = ref_dir, 
        scaffold = "data/assemblies/{sample}/scaffolds.fasta" 
    output: 
        "data/scaffolds/{sample}_scaffolds.fasta" 
    log: 
        "data/scaffolds/logs/{sample}.log"
    message: "Executing assembly refinement using Medusa."
    shell:
        "cp -r medusa/medusa_scripts . && java -jar medusa/medusa.jar -f {input.ref} -i {input.scaffold} -o {output} -v > {log}"

