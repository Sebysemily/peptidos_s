configfile: "config/config.yml"

CURATED_FASTAS = config["curated_fastas"]
PEPTIDE_SETS = list(CURATED_FASTAS.keys())
BATCHING = config.get("batching", {})


def n_batches(peptide_set):
    batches = int(BATCHING.get(peptide_set, 1))
    if batches < 1:
        raise ValueError(f"batching for {peptide_set} must be at least 1")
    return batches


def batch_ids_for(peptide_set):
    return [str(i) for i in range(1, n_batches(peptide_set) + 1)]


include: "rules/pre_processing.smk"
include: "rules/tox_check.smk"

PRE_PROCESSING_TARGETS = [
    *expand(
        "data/curated_md-lais/mmseqs2/{peptide_set}/clusters_{peptide_set}_rep_seq.fasta",
        peptide_set=PEPTIDE_SETS,
    )
    + expand(
        "data/curated_md-lais/mmseqs2/{peptide_set}/clusters_{peptide_set}_cluster.tsv",
        peptide_set=PEPTIDE_SETS,
    )
]

TOX_CHECK_TARGETS = [
    *expand(
        "results/toxinpred3/{peptide_set}/clusters_{peptide_set}_rep_seq_toxinpred3.csv",
        peptide_set=PEPTIDE_SETS,
    )
]

TOXTELLER_TARGETS = [
    *expand(
        "results/toxteller/{peptide_set}/clusters_{peptide_set}_rep_seq_toxteller.csv",
        peptide_set=PEPTIDE_SETS,
    )
]


rule all:
    input:
        PRE_PROCESSING_TARGETS + TOX_CHECK_TARGETS + TOXTELLER_TARGETS
