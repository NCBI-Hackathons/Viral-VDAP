import Bio.AlignIO

##### Helper Functions #####
## Modified from Broad Institutes viral-ngs pipeline ##

def fastaMaker(seq, idVal, linewidth=60):
    """Generate FASTA file""" 
    assert linewidth > 0
        
    yield ">{}\n".format(idVal)

    while len(seq) > linewidth:
        line = seq[:linewidth]
        seq = seq[linewidth:]
        yield "{}\n".format(line.upper())

    if seq:
        yield seq.upper() + "\n"

def makeFastaFile(seqs, idVal, outFasta):
    """Write to FASTA file"""
   
    with open(outFasta, 'wt') as outf:
        for line in fastaMaker(seqs, idVal):
            outf.write(line)

    return outFasta

def modify_scaffold(aln_file, idVal, outFasta):
    """Modifies input scaffold by trimming to the length of the reference"""

    # read in pairwise alignment file
    aln = Bio.AlignIO.read(str(aln_file), "fasta")
    
    if len(aln) !=2: 
        raise Exception("Alignment does not contain exactly 2 sequences, %s found" % len(aln))
    else:
        ref_idx = 0 # reference is first sequence
        consensus_idx = 1 # sample is second sequence
    
    # get sequences
    ref = str(aln[ref_idx].seq)
    consensus = str(aln[consensus_idx].seq)
       
    if len(ref) != len(consensus):
        raise Exception("Improper alignment")
    else:
        # convert to list
        ref = list(ref)
        consensus = list(consensus)
 
    # trim ends of consensus so it doesn't go beyond given reference genome
    for end_iterator in (range(len(ref)), reversed(range(len(ref)))):
        for i in end_iterator:
            if ref[i] != "-":
                break
            else:
                consensus[i] = "-"

    # fill out ends of the consensus with reference sequence
    for end_iterator in (range(len(ref)), reversed(range(len(ref)))):
        for i in end_iterator:
            if consensus[i] != "-":
                break
            else:
                consensus[i] = ref[i]

    # remove gaps within consensus
    final_consensus = "".join(consensus).replace("-","")

    # write to file
    return makeFastaFile(final_consensus, idVal, str(outFasta)) 


##### Rule #####
rule merge_contigs:
    """Merge contigs in sample's scaffolds FASTA file and concatenate with reference FASTA file"""
    input:
        ref = ref_fasta,
        scaffolds = "data/final_scaffolds/{sample}_scaffolds.fasta"
    output:
        temp_file = "data/final_scaffolds/{sample}_merged.fasta",
        final_file = "data/final_scaffolds/ref_and_{sample}.fasta"
    params: "{sample}"        
    message: "Preparing FASTA file for pairwise alignment."
    shell:
        """grep -v "^>" {input.scaffolds} > {output.temp_file} && echo -e ">{params}\n$(cat {output.temp_file})" > {output.temp_file} && cat {input.ref} {output.temp_file} > {output.final_file}"""

 
##### Rule #####
rule mafft:
    """Pairwise alignment of reference and sample sequences with MAFFT"""
    input:
        "data/final_scaffolds/ref_and_{sample}.fasta"
    output:
        "data/pairwise_alignments/{sample}_pairwise_aln.fasta"
    params:
        "--auto --globalpair --reorder --maxiterate 1000"
    message: "Executing pairwise alignment with MAFFT."
    log: 
        "data/pairwise_alignments/logs/{sample}.log"
    threads: num_threads
    shell: 
        "(mafft --thread {threads} {params} {input} > {output}) 2> {log}"

##### Rule #####
rule create_consensus:
    """Modify sample sequence in pairwise alignment file to create consensus sequence FASTA file"""
    input:
        "data/pairwise_alignments/{sample}_pairwise_aln.fasta"
    output:
        "data/consensus_fasta/{sample}/{sample}_consensus.fasta"
    params: "{sample}"
    message: "Creating consensus FASTA file."
    run: 
        modify_scaffold(input,params,output) 
