configfile: "config/config.yml"

CURATED_FASTAS = config["curated_fastas"]
PEPTIDE_SETS = list(CURATED_FASTAS.keys())
BATCHING = config.get("batching", {})


def n_batches(peptide_set):
    batches = int(BATCHING.get(peptide_set, 1))
    if batches < 1:
        raise ValueError(f"batching for {peptide_set} must be at least 1")
    return batches


include: "rules/pre_processing.smk"
include: "rules/tox_check.smk"
include: "rules/hem_check.smk"

SETUP_TARGETS = [
    "results/setup/.external_resources_checked",
]

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
        "results/tox_check/toxinpred3/{peptide_set}/clusters_{peptide_set}_rep_seq_toxinpred3.csv",
        peptide_set=PEPTIDE_SETS,
    )
]

TOXTELLER_TARGETS = [
    *expand(
        "results/tox_check/toxteller/{peptide_set}/clusters_{peptide_set}_rep_seq_toxteller.csv",
        peptide_set=PEPTIDE_SETS,
    )
]

CAPTP_TARGETS = [
    *expand(
        "results/tox_check/captp/{peptide_set}/clusters_{peptide_set}_rep_seq_captp.csv",
        peptide_set=PEPTIDE_SETS,
    )
]

TOXICITY_SUMMARY_TARGETS = [
    *expand(
        "results/tox_check/toxicity_summary/{peptide_set}/clusters_{peptide_set}_toxicity_summary.csv",
        peptide_set=PEPTIDE_SETS,
    )
]

HEMOPI2_TARGETS = [
    *expand(
        (
            "results/hemo_check/hemopi2_classification/{peptide_set}/"
            "clusters_{peptide_set}_rep_seq_hemopi2_classification.csv"
        ),
        peptide_set=PEPTIDE_SETS,
    ),
    *expand(
        (
            "results/hemo_check/hemopi2_regression/{peptide_set}/"
            "clusters_{peptide_set}_rep_seq_hemopi2_regression.csv"
        ),
        peptide_set=PEPTIDE_SETS,
    )
]


rule all:
    input:
        (
            SETUP_TARGETS
            + PRE_PROCESSING_TARGETS
            + TOX_CHECK_TARGETS
            + TOXTELLER_TARGETS
            + CAPTP_TARGETS
            + TOXICITY_SUMMARY_TARGETS
            + HEMOPI2_TARGETS
        )
