# path to reference directory
ref_dir = config["ref"]["dir"]

##### Rule #####
rule refine_assembly:
    """Refine assembly using Medusa"""
    input: 
        ref = ref_dir, 
        scaffold = "data/assemblies/{sample}/scaffolds.fasta" 
    output: 
        "data/final_scaffolds/{sample}_scaffolds.fasta" 
    log: 
        "data/final_scaffolds/logs/{sample}.log"
    message: "Executing assembly refinement using Medusa."
    shell:
        "cp -r medusa/medusa_scripts . && java -jar medusa/medusa.jar -f {input.ref} -i {input.scaffold} -o {output} -v > {log}"

