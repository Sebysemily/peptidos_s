import re

THREADS = int(config.get("max_threads", 8))

MMSEQS_CONFIG = config.get("mmseqs", {})
MMSEQS_MIN_SEQ_ID = MMSEQS_CONFIG.get("min_seq_id", 0.95)
MMSEQS_COVERAGE = MMSEQS_CONFIG.get("coverage", 0.8)
MMSEQS_COV_MODE = int(MMSEQS_CONFIG.get("cov_mode", 1))
MMSEQS_SENSITIVITY = MMSEQS_CONFIG.get("sensitivity", 7.5)
MMSEQS_SUB_MAT = MMSEQS_CONFIG.get(
    "sub_mat",
    "aa:blosum62.out,nucl:nucleotide.out",
)

TOXTELLER_PROGRAM_DIR = "resources/ToxTeller"
CAPTP_PROGRAM_DIR = "resources/CAPTP"

PEPTIDE_FASTAS = config["curated_fastas"]
PEPTIDE_SET = "|".join(re.escape(name) for name in PEPTIDE_FASTAS)
wildcard_constraints:
    peptide_set=PEPTIDE_SET


rule check_external_resources:
    output:
        touch("results/setup/.external_resources_checked")
    priority: 100
    params:
        toxteller=TOXTELLER_PROGRAM_DIR,
        captp=CAPTP_PROGRAM_DIR,
    shell:
        r"""
        missing=0

        if [ ! -f "{params.toxteller}/program_resource/toxteller.py" ]; then
            echo "ToxTeller checkout is missing or incomplete: {params.toxteller}" >&2
            echo "Expected file: {params.toxteller}/program_resource/toxteller.py" >&2
            missing=1
        fi

        if [ ! -f "{params.captp}/main.py" ]; then
            echo "CAPTP checkout is missing or incomplete: {params.captp}" >&2
            echo "Expected file: {params.captp}/main.py" >&2
            missing=1
        fi

        if [ "$missing" -ne 0 ]; then
            echo "Initialize external resources with:" >&2
            echo "  git submodule update --init --recursive" >&2
            exit 1
        fi
        """


rule mmseqs_easy_cluster:
    input:
        fasta=lambda wildcards: PEPTIDE_FASTAS[wildcards.peptide_set],
    output:
        rep_seq="data/curated_md-lais/mmseqs2/{peptide_set}/clusters_{peptide_set}_rep_seq.fasta",
        cluster_tsv="data/curated_md-lais/mmseqs2/{peptide_set}/clusters_{peptide_set}_cluster.tsv",
    params:
        prefix="data/curated_md-lais/mmseqs2/{peptide_set}/clusters_{peptide_set}",
        tmp="data/curated_md-lais/mmseqs2/{peptide_set}/tmp_{peptide_set}",
        outdir="data/curated_md-lais/mmseqs2/{peptide_set}",
        min_seq_id=MMSEQS_MIN_SEQ_ID,
        coverage=MMSEQS_COVERAGE,
        cov_mode=MMSEQS_COV_MODE,
        sensitivity=MMSEQS_SENSITIVITY,
        sub_mat=MMSEQS_SUB_MAT,
    threads: THREADS
    conda:
        "../envs/pre_processing.yml"
    shell:
        r"""
        mkdir -p {params.outdir} {params.tmp}
        mmseqs easy-cluster {input.fasta} {params.prefix} {params.tmp} \
            --min-seq-id {params.min_seq_id} \
            -c {params.coverage} \
            --cov-mode {params.cov_mode} \
            -s {params.sensitivity} \
            --sub-mat {params.sub_mat} \
            --threads {threads}
        """
