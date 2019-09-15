import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# Config file:
configfile: "config.yaml"

# Read in samples.tsv file
samples = pd.read_table(config["samples"]).set_index("sample",drop=False)

# Path to adapters file
adapters = config["adapters"]

# Number of threads to use
num_threads = int(config["threads"])

# Check for platform: iontorrent or illumina
platform = config["spades"]["platform"]

# Get list of k-mer sizes
kmers = config["spades"]["kmer_sizes"]

# if platform is iontorrent, create parameter to be added to SPAdes
if platform == "iontorrent":
    is_ion = "--iontorrent"
else:
    is_ion = ""

# if kmer_sizes is not empty, create parameter to be added to SPAdes
if kmers is not None:
    has_kmers = "-k {}".format(kmers)
else:
    has_kmers = ""

# combine parameters
add_ons = is_ion + " " + has_kmers

# path to reference directory
ref_dir = config["ref"]["dir"]

# path to reference genome
ref = config["ref"]["genome"]

##### Helper functions #####

def get_fastq(wildcards):
    """Get FASTQ files of given sample."""
    return samples.loc[wildcards.sample, ["fq1", "fq2"]].dropna()

def get_id(fasta_file):
    """Get sequence id from FASTA file for longest sequence"""

    max_len = 0
    max_id = ""

    for record in SeqIO.parse(open(str(fasta_file),mode="r"),"fasta"):
      if len(record) > max_len:
        max_len = len(record)
        max_id = record.name

    return max_id
 
