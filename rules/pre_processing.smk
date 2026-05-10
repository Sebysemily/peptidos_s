MMSEQS_CONFIG = config.get("mmseqs", {})
MMSEQS_THREADS = int(MMSEQS_CONFIG.get("threads", config.get("max_threads", 8)))
MMSEQS_MIN_SEQ_ID = MMSEQS_CONFIG.get("min_seq_id", 0.95)
MMSEQS_COVERAGE = MMSEQS_CONFIG.get("coverage", 0.8)
MMSEQS_COV_MODE = int(MMSEQS_CONFIG.get("cov_mode", 1))
MMSEQS_SENSITIVITY = MMSEQS_CONFIG.get("sensitivity", 7.5)
MMSEQS_SUB_MAT = MMSEQS_CONFIG.get(
    "sub_mat",
    "aa:blosum62.out,nucl:nucleotide.out",
)

PEPTIDE_FASTAS = {
    "10_25": "data/curated_md-lais/10_25/alls_RangeSelected_10_25_Unique.fasta",
    "25_50": "data/curated_md-lais/25_50/alls_RangeSelected_25_50_curado.fasta",
}


wildcard_constraints:
    peptide_range="10_25|25_50"


rule mmseqs_easy_cluster:
    input:
        fasta=lambda wildcards: PEPTIDE_FASTAS[wildcards.peptide_range],
    output:
        rep_seq="data/curated_md-lais/mmseqs2/{peptide_range}/clusters_{peptide_range}_rep_seq.fasta",
        cluster_tsv="data/curated_md-lais/mmseqs2/{peptide_range}/clusters_{peptide_range}_cluster.tsv",
    params:
        prefix="data/curated_md-lais/mmseqs2/{peptide_range}/clusters_{peptide_range}",
        tmp="data/curated_md-lais/mmseqs2/{peptide_range}/tmp_{peptide_range}",
        outdir="data/curated_md-lais/mmseqs2/{peptide_range}",
        min_seq_id=MMSEQS_MIN_SEQ_ID,
        coverage=MMSEQS_COVERAGE,
        cov_mode=MMSEQS_COV_MODE,
        sensitivity=MMSEQS_SENSITIVITY,
        sub_mat=MMSEQS_SUB_MAT,
    threads: MMSEQS_THREADS
    conda:
        "../envs/mmseq2.yml"
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
