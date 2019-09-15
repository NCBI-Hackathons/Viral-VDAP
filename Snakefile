include: "rules/common.smk"

##### Target rules #####
rule all:
    input:
        expand("data/vcf/{sample}/{sample}.vcf", sample=samples["sample"])


##### Modules #####

include: "rules/trim.smk"
include: "rules/denovo_assembly.smk"
include: "rules/refine_assembly.smk"
include: "rules/qc.smk"
include: "rules/consensus.smk"
include: "rules/map_reads.smk"
include: "rules/variant_calling.smk"
