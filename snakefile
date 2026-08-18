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
NON_HEMO_FASTA_ROOT = config.get(
    "non_hemo_fasta_root", "data/derived/non_hemo"
)
POST_CHECKS_FASTA_ROOT = config.get(
    "post_checks_fasta_root", "data/derived/post_checks"
)


def n_batches(peptide_set):
    batches = int(BATCHING.get(peptide_set, 1))
    if batches < 1:
        raise ValueError(f"batching for {peptide_set} must be at least 1")
    return batches


include: "rules/pre_processing.smk"
include: "rules/tox_check.smk"
include: "rules/hem_check.smk"
include: "rules/inmuno_check.smk"
include: "rules/properties.smk"
include: "rules/ACP_predictors.smk" # ACPScanner disabled internally, but file has AntiCP2

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

# AMPTOX_TARGETS = [
#     *expand(
#         "results/tox_check/amptox/{peptide_set}/clusters_{peptide_set}_rep_seq_amptox.csv",
#         peptide_set=PEPTIDE_SETS,
#     )
# ]

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


HEMO_DL_TARGETS = [
    *expand(
        (
            "results/hemo_check/hemo_DL/{peptide_set}/"
            "clusters_{peptide_set}_rep_seq_hemo_DL.csv"
        ),
        peptide_set=PEPTIDE_SETS,
    )
]

HEMOLYTIC_SUMMARY_TARGETS = [
    *expand(
        (
            "results/hemo_check/hemolytic_summary/{peptide_set}/"
            "clusters_{peptide_set}_hemolytic_summary.csv"
        ),
        peptide_set=PEPTIDE_SETS,
    )
]

NON_HEMO_FASTA_TARGETS = [
    *expand(
        (
            f"{NON_HEMO_FASTA_ROOT}/{{peptide_set}}/"
            "clusters_{peptide_set}_rep_seq_non_hemo.fasta"
        ),
        peptide_set=PEPTIDE_SETS,
    ),
    *expand(
        f"{NON_HEMO_FASTA_ROOT}/{{peptide_set}}/batches",
        peptide_set=PEPTIDE_SETS,
    ),
]

ALGPRED2_TARGETS = [
    *expand(
        (
            "results/inmuno_check/algpred2/{peptide_set}/"
            "clusters_{peptide_set}_rep_seq_algpred2.csv"
        ),
        peptide_set=PEPTIDE_SETS,
    )
]

ALLERGENAI_TARGETS = [
    *expand(
        (
            "results/inmuno_check/allergenai/{peptide_set}/"
            "clusters_{peptide_set}_rep_seq_allergenai.csv"
        ),
        peptide_set=[p for p in PEPTIDE_SETS if p != "10_25"],
    )
]

ALLERTRANS_TARGETS = [
    *expand(
        (
            "results/inmuno_check/allertrans/{peptide_set}/"
            "clusters_{peptide_set}_rep_seq_allertrans.csv"
        ),
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


INMUNO_SUMMARY_TARGETS = [
    *expand(
        "results/inmuno_check/inmuno_summary/{peptide_set}/clusters_{peptide_set}_inmuno_summary.csv",
        peptide_set=PEPTIDE_SETS,
    )
]

POST_CHECKS_FASTA_TARGETS = [
    *expand(
        f"{POST_CHECKS_FASTA_ROOT}/{{peptide_set}}/clusters_{{peptide_set}}_rep_seq_post_checks.fasta",
        peptide_set=PEPTIDE_SETS,
    )
]

ANTICP2_TARGETS = [
    *expand(
        "results/acp_predictors/anticp2/{peptide_set}/clusters_{peptide_set}_rep_seq_anticp2.csv",
        peptide_set=PEPTIDE_SETS,
    )
]

# ACP_PREDICTORS_TARGETS = [
#     *expand(
#         "results/acp_predictors/acpscanner/{peptide_set}/clusters_{peptide_set}_rep_seq_acpscanner.csv",
#         peptide_set=PEPTIDE_SETS,
#     )
# ]

ACP_OPE_TARGETS = [
    *expand(
        "results/acp_predictors/acp_ope/{peptide_set}/clusters_{peptide_set}_rep_seq_acp_ope.csv",
        peptide_set=PEPTIDE_SETS,
    )
]

PROCESS_GROUPS_TOGETHER = config.get("properties", {}).get("process_groups_together", True)
if PROCESS_GROUPS_TOGETHER:
    PROPERTIES_TARGETS = ["metadata/characteristics.csv", "metadata/filtering.csv"]
else:
    PROPERTIES_TARGETS = [
        *expand(
            "metadata/{peptide_set}/characteristics.csv",
            peptide_set=PEPTIDE_SETS,
        ),
        "metadata/filtering.csv"
    ]

rule all:
    input:
        (
            SETUP_TARGETS
            + PRE_PROCESSING_TARGETS
            + TOX_CHECK_TARGETS
            + TOXTELLER_TARGETS
            + CAPTP_TARGETS
            # + AMPTOX_TARGETS
            + TOXICITY_SUMMARY_TARGETS
            + NON_TOXIC_FASTA_TARGETS
            + HEMOPI2_TARGETS
            + MACREL_TARGETS
            + HEPAD_TARGETS
            + HEMO_DL_TARGETS
            + HEMOLYTIC_SUMMARY_TARGETS
            + NON_HEMO_FASTA_TARGETS
            + ALGPRED2_TARGETS
            + ALLERGENAI_TARGETS
            + ALLERTRANS_TARGETS
            + INMUNO_SUMMARY_TARGETS
            + POST_CHECKS_FASTA_TARGETS
            + PROPERTIES_TARGETS
            + ANTICP2_TARGETS
            + ACP_OPE_TARGETS
            # + ACP_PREDICTORS_TARGETS # Disabled: requires SPIDER3, ESM, and PDB structures
        )


# HEPAD-only entry point: skips tox re-runs when non-toxic batches already exist.
# Usage:
#   snakemake hepad_all --rerun-triggers mtime -j <N> --use-conda
rule hepad_all:
    input:
        HEPAD_ONLY_TARGETS


