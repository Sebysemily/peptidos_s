configfile: "config/config.yml"

import importlib.util
# Load the helper module from code/ directly using importlib to avoid conflict with standard library 'code' module
spec = importlib.util.spec_from_file_location("paralelism_based_on_total_ram", "code/paralelism_based_on_total_ram.py")
helper = importlib.util.module_from_spec(spec)
spec.loader.exec_module(helper)
get_mem_mb = helper.get_mem_mb
get_hepad_mem_mb = helper.get_hepad_mem_mb
set_workflow_resources = helper.set_workflow_resources

set_workflow_resources(workflow, leeway_gb=10)





CURATED_FASTAS = config["curated_fastas"]
PEPTIDE_SETS = list(CURATED_FASTAS.keys())
BATCHING = config.get("batching", {})
BATCH_FASTA_ROOT = config.get("batch_fasta_root", "data/derived/batches")
TOX_CHECK_BATCH_DIR = f"{BATCH_FASTA_ROOT}/tox_check"
TOXTELLER_BATCH_DIR = f"{BATCH_FASTA_ROOT}/toxteller"
CAPTP_BATCH_DIR = f"{BATCH_FASTA_ROOT}/captp"
NON_TOXIC_FASTA_ROOT = config.get(
    "non_toxic_fasta_root", "data/derived/non_toxic"
)


def n_batches(peptide_set):
    batches = int(BATCHING.get(peptide_set, 1))
    if batches < 1:
        raise ValueError(f"batching for {peptide_set} must be at least 1")
    return batches


include: "rules/pre_processing.smk"
include: "rules/tox_check.smk"
include: "rules/hem_check.smk"

SETUP_TARGETS = [
    ".snakemake/checks/external_resources_checked",
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
    ),
    *expand(
        "results/tox_check/toxinpred3/{peptide_set}/.batches_validated",
        peptide_set=PEPTIDE_SETS,
    )
]

TOXTELLER_TARGETS = [
    *expand(
        "results/tox_check/toxteller/{peptide_set}/clusters_{peptide_set}_rep_seq_toxteller.csv",
        peptide_set=PEPTIDE_SETS,
    ),
    *expand(
        "results/tox_check/toxteller/{peptide_set}/.batches_validated",
        peptide_set=PEPTIDE_SETS,
    )
]

CAPTP_TARGETS = [
    *expand(
        "results/tox_check/captp/{peptide_set}/clusters_{peptide_set}_rep_seq_captp.csv",
        peptide_set=PEPTIDE_SETS,
    ),
    *expand(
        "results/tox_check/captp/{peptide_set}/.batches_validated",
        peptide_set=PEPTIDE_SETS,
    )
]

TOXICITY_SUMMARY_TARGETS = [
    *expand(
        "results/tox_check/toxicity_summary/{peptide_set}/clusters_{peptide_set}_toxicity_summary.csv",
        peptide_set=PEPTIDE_SETS,
    )
]

NON_TOXIC_FASTA_TARGETS = [
    *expand(
        (
            f"{NON_TOXIC_FASTA_ROOT}/{{peptide_set}}/"
            "clusters_{peptide_set}_rep_seq_non_toxic.fasta"
        ),
        peptide_set=PEPTIDE_SETS,
    ),
    *expand(
        f"{NON_TOXIC_FASTA_ROOT}/{{peptide_set}}/batches",
        peptide_set=PEPTIDE_SETS,
    ),
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

MACREL_TARGETS = [
    *expand(
        (
            "results/hemo_check/macrel/{peptide_set}/"
            "clusters_{peptide_set}_rep_seq_macrel.csv"
        ),
        peptide_set=PEPTIDE_SETS,
    )
]

HEPAD_DATASETS = config.get(
    "hepad_datasets", ["Hmp1", "Hmp2", "Hmp3", "Hmpm"]
)
HEPAD_TARGETS = [
    *expand(
        (
            "results/hemo_check/hepad/{hepad_dataset}/{peptide_set}/"
            "clusters_{peptide_set}_rep_seq_hepad_{hepad_dataset}.csv"
        ),
        hepad_dataset=HEPAD_DATASETS,
        peptide_set=PEPTIDE_SETS,
    )
]


HEPAD_ONLY_TARGETS = [
    *expand(
        f"{NON_TOXIC_FASTA_ROOT}/{{peptide_set}}/batches",
        peptide_set=PEPTIDE_SETS,
    ),
    *HEPAD_TARGETS,
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
            + NON_TOXIC_FASTA_TARGETS
            + HEMOPI2_TARGETS
            + MACREL_TARGETS
            + HEPAD_TARGETS
        )


# HEPAD-only entry point: skips tox re-runs when non-toxic batches already exist.
# Usage:
#   snakemake hepad_all --rerun-triggers mtime -j <N> --use-conda
rule hepad_all:
    input:
        HEPAD_ONLY_TARGETS
