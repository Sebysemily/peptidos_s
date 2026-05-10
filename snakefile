configfile: "config/config.yml"

include: "rules/pre_processing.smk"

PEPTIDE_RANGES = ["10_25", "25_50"]
PRE_PROCESSING_TARGETS = [
    *expand(
        "data/curated_md-lais/mmseqs2/{peptide_range}/clusters_{peptide_range}_rep_seq.fasta",
        peptide_range=PEPTIDE_RANGES,
    )
    + expand(
        "data/curated_md-lais/mmseqs2/{peptide_range}/clusters_{peptide_range}_cluster.tsv",
        peptide_range=PEPTIDE_RANGES,
    )
]


rule all:
    input:
        PRE_PROCESSING_TARGETS
