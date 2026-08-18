rule download_acpscanner_resources:
    output:
        "resources/acpscanner/requirements.txt"
    shell:
        r"""
        mkdir -p resources/acpscanner
        wget -q http://acpscanner.denglab.org/test/materials.zip -O resources/acpscanner/materials.zip
        cd resources/acpscanner
        unzip -q -o materials.zip
        # Move files from material folder to resources/acpscanner
        mv material/* .
        rmdir material
        rm materials.zip
        """

rule run_acpscanner:
    input:
        fasta=f"{POST_CHECKS_FASTA_ROOT}/{{peptide_set}}/clusters_{{peptide_set}}_rep_seq_post_checks.fasta",
        resources="resources/acpscanner/requirements.txt"
    output:
        csv="results/acp_predictors/acpscanner/{peptide_set}/clusters_{peptide_set}_rep_seq_acpscanner.csv"
    conda:
        "../envs/acp_predictors.yml"
    shell:
        r"""
        python code/acp_predictors/run_acpscanner.py \
            --fasta {input.fasta} \
            --resources resources/acpscanner \
            --output {output.csv}
        """

checkpoint split_post_checks_batches:
    input:
        script=ancient("code/split_fasta_batches.py"),
        fasta=f"{POST_CHECKS_FASTA_ROOT}/{{peptide_set}}/clusters_{{peptide_set}}_rep_seq_post_checks.fasta",
    output:
        directory(f"{POST_CHECKS_FASTA_ROOT}/{{peptide_set}}/batches"),
    params:
        n_batches=lambda wildcards: n_batches(wildcards.peptide_set),
        id_prefix=lambda wildcards: f"post_checks_{wildcards.peptide_set}",
    conda:
        "../envs/tox_check/toxinpred3_captp.yml"
    shell:
        r"""
        python {input.script} \
            --input {input.fasta} \
            --outdir {output} \
            --n-batches {params.n_batches} \
            --id-prefix {params.id_prefix}
        """

def post_checks_batch_dir(wildcards):
    return checkpoints.split_post_checks_batches.get(
        peptide_set=wildcards.peptide_set
    ).output[0]

def post_checks_batch_fasta(wildcards):
    return ancient(
        f"{post_checks_batch_dir(wildcards)}/batch_{wildcards.batch_id}.fasta"
    )

def post_checks_batch_mapping(wildcards):
    return ancient(
        f"{post_checks_batch_dir(wildcards)}/batch_{wildcards.batch_id}.mapping.csv"
    )

def post_checks_batch_reports(wildcards, tool):
    batch_dir = post_checks_batch_dir(wildcards)
    batch_ids = glob_wildcards(
        f"{batch_dir}/batch_{{batch_id}}.fasta"
    ).batch_id
    batch_ids = sorted(batch_ids, key=int)
    reports = expand(
        "results/acp_predictors/{tool}/{peptide_set}/batches/batch_{batch_id}.csv",
        tool=tool,
        peptide_set=wildcards.peptide_set,
        batch_id=batch_ids,
    )
    return [ancient(report) for report in reports]

rule anticp2_batch:
    input:
        fasta=post_checks_batch_fasta,
        mapping=post_checks_batch_mapping,
    output:
        report=(
            "results/acp_predictors/anticp2/{peptide_set}/batches/"
            "batch_{batch_id}.csv"
        ),
    params:
        outdir="results/acp_predictors/anticp2/{peptide_set}/batches",
        raw_report="results/acp_predictors/anticp2/{peptide_set}/batches/raw_batch_{batch_id}.txt",
    threads: 1
    resources:
        mem_mb=get_mem_mb
    conda:
        "../envs/anticp2.yml"
    shell:
        r"""
        mkdir -p {params.outdir}
        python code/acp_predictors/run_anticp2.py \
            --fasta {input.fasta} \
            --mapping {input.mapping} \
            --raw-output {params.raw_report} \
            --csv-output {output.report}
        """

rule merge_anticp2_batches:
    input:
        lambda wildcards: post_checks_batch_reports(wildcards, "anticp2")
    output:
        report=(
            "results/acp_predictors/anticp2/{peptide_set}/"
            "clusters_{peptide_set}_rep_seq_anticp2.csv"
        ),
        validated=touch(
            "results/acp_predictors/anticp2/{peptide_set}/.batches_validated"
        ),
    conda:
        "../envs/anticp2.yml"
    shell:
        r"""
        python code/merge_csv_reports.py \
            --inputs {input} \
            --output {output.report}
        """

rule acp_ope_batch:
    input:
        fasta=post_checks_batch_fasta,
        mapping=post_checks_batch_mapping,
    output:
        report=(
            "results/acp_predictors/acp_ope/{peptide_set}/batches/"
            "batch_{batch_id}.csv"
        ),
    params:
        outdir="results/acp_predictors/acp_ope/{peptide_set}/batches",
        model_dir="resources/acp-ope/models",
    threads: 1
    resources:
        mem_mb=get_mem_mb
    conda:
        "../envs/acp_ope.yml"
    shell:
        r"""
        mkdir -p {params.outdir}
        python code/acp_predictors/run_acp_ope.py \
            --fasta {input.fasta} \
            --mapping {input.mapping} \
            --model-dir {params.model_dir} \
            --csv-output {output.report}
        """

rule merge_acp_ope_batches:
    input:
        lambda wildcards: post_checks_batch_reports(wildcards, "acp_ope")
    output:
        report=(
            "results/acp_predictors/acp_ope/{peptide_set}/"
            "clusters_{peptide_set}_rep_seq_acp_ope.csv"
        ),
        validated=touch(
            "results/acp_predictors/acp_ope/{peptide_set}/.batches_validated"
        ),
    conda:
        "../envs/acp_ope.yml"
    shell:
        r"""
        python code/merge_csv_reports.py \
            --inputs {input} \
            --output {output.report}
        """
