import re

THREADS = int(config.get("max_threads", 8))

PEPTIDE_SET = "|".join(
    re.escape(name) for name in config["curated_fastas"]
)


wildcard_constraints:
    peptide_set=PEPTIDE_SET,
    batch_id=r"\d+"


checkpoint split_tox_check_batches:
    input:
        script=ancient("code/split_fasta_batches.py"),
        rep_seq=(
            "data/curated_md-lais/mmseqs2/{peptide_set}/"
            "clusters_{peptide_set}_rep_seq.fasta"
        ),
    output:
        directory(f"{TOX_CHECK_BATCH_DIR}/{{peptide_set}}"),
    params:
        n_batches=lambda wildcards: n_batches(wildcards.peptide_set),
        id_prefix=lambda wildcards: f"tox_check_{wildcards.peptide_set}",
    conda:
        "../envs/tox_check/toxinpred3_captp.yml"
    shell:
        r"""
        python {input.script} \
            --input {input.rep_seq} \
            --outdir {output} \
            --n-batches {params.n_batches} \
            --id-prefix {params.id_prefix}
        """


def tox_check_batch_dir(wildcards):
    return checkpoints.split_tox_check_batches.get(
        peptide_set=wildcards.peptide_set
    ).output[0]


def tox_check_batch_fasta(wildcards):
    return ancient(
        f"{tox_check_batch_dir(wildcards)}/batch_{wildcards.batch_id}.fasta"
    )


def tox_check_batch_mapping(wildcards):
    return ancient(
        f"{tox_check_batch_dir(wildcards)}/batch_{wildcards.batch_id}.mapping.csv"
    )


def tox_check_batch_reports(wildcards, tool):
    batch_dir = tox_check_batch_dir(wildcards)
    batch_ids = glob_wildcards(
        f"{batch_dir}/batch_{{batch_id}}.fasta"
    ).batch_id
    batch_ids = sorted(batch_ids, key=int)
    reports = expand(
        "results/tox_check/{tool}/{peptide_set}/batches/batch_{batch_id}.csv",
        tool=tool,
        peptide_set=wildcards.peptide_set,
        batch_id=batch_ids,
    )
    return [ancient(report) for report in reports]


checkpoint split_toxteller_batches:
    input:
        script=ancient("code/split_fasta_batches.py"),
        rep_seq=(
            "data/curated_md-lais/mmseqs2/{peptide_set}/"
            "clusters_{peptide_set}_rep_seq.fasta"
        ),
    output:
        directory(f"{TOXTELLER_BATCH_DIR}/{{peptide_set}}"),
    params:
        max_sequences_per_batch=9500,
    conda:
        "../envs/tox_check/toxteller.yml"
    shell:
        r"""
        python {input.script} \
            --input {input.rep_seq} \
            --outdir {output} \
            --max-records-per-batch {params.max_sequences_per_batch}
        """


def toxteller_batch_dir(wildcards):
    return checkpoints.split_toxteller_batches.get(
        peptide_set=wildcards.peptide_set
    ).output[0]


def toxteller_batch_fasta(wildcards):
    return ancient(
        f"{toxteller_batch_dir(wildcards)}/batch_{wildcards.batch_id}.fasta"
    )


def toxteller_batch_reports(wildcards):
    batch_dir = toxteller_batch_dir(wildcards)
    batch_ids = glob_wildcards(
        f"{batch_dir}/batch_{{batch_id}}.fasta"
    ).batch_id
    batch_ids = sorted(batch_ids, key=int)
    reports = expand(
        "results/tox_check/toxteller/{peptide_set}/batches/batch_{batch_id}.csv",
        peptide_set=wildcards.peptide_set,
        batch_id=batch_ids,
    )
    return [ancient(report) for report in reports]


rule toxinpred3_batch:
    input:
        fasta=tox_check_batch_fasta,
        mapping=tox_check_batch_mapping,
    output:
        report=(
            "results/tox_check/toxinpred3/{peptide_set}/batches/"
            "batch_{batch_id}.csv"
        ),
    params:
        outdir="results/tox_check/toxinpred3/{peptide_set}/batches",
        workdir="results/tox_check/toxinpred3/{peptide_set}/work/batch_{batch_id}",
        raw_report="results/tox_check/toxinpred3/{peptide_set}/work/batch_{batch_id}/raw_toxinpred3.csv",
    threads: 1
    resources:
        mem_mb=get_mem_mb
    conda:
        "../envs/tox_check/toxinpred3_captp.yml"
    shell:
        r"""
        mkdir -p {params.outdir}
        rm -rf {params.workdir}
        mkdir -p {params.workdir}
        input_fasta="$(realpath {input.fasta})"
        raw_report="$(realpath -m {params.raw_report})"
        output_report="$(realpath -m {output.report})"
        cd {params.workdir}
        toxinpred3 -i "$input_fasta" -o "$raw_report" -m 2 -d 2
        cd -
        python code/tox_check/annotate_toxinpred3_report.py \
            --raw-report "$raw_report" \
            --mapping {input.mapping} \
            --output "$output_report"
        """


rule merge_toxinpred3_batches:
    input:
        lambda wildcards: tox_check_batch_reports(wildcards, "toxinpred3")
    output:
        report=(
            "results/tox_check/toxinpred3/{peptide_set}/"
            "clusters_{peptide_set}_rep_seq_toxinpred3.csv"
        ),
        validated=touch(
            "results/tox_check/toxinpred3/{peptide_set}/.batches_validated"
        ),
    conda:
        "../envs/tox_check/toxinpred3_captp.yml"
    shell:
        r"""
        python code/merge_csv_reports.py \
            --inputs {input} \
            --output {output.report}
        """

rule toxteller_batch:
    input:
        resources_checked=ancient(".snakemake/checks/external_resources_checked"),
        fasta=toxteller_batch_fasta,
    output:
        report=(
            "results/tox_check/toxteller/{peptide_set}/batches/"
            "batch_{batch_id}.csv"
        ),
    params:
        outdir="results/tox_check/toxteller/{peptide_set}/batches",
        tool_dir=TOXTELLER_PROGRAM_DIR,
    threads: 1
    resources:
        mem_mb=get_mem_mb
    conda:
        "../envs/tox_check/toxteller.yml"
    shell:
        r"""
        mkdir -p {params.outdir}
        input_fasta="$(realpath {input.fasta})"
        output_report="$(realpath -m {output.report})"
        cd {params.tool_dir}/program_resource
        python toxteller.py "$input_fasta"
        mv "${{input_fasta}}.csv" "$output_report"
        """


rule merge_toxteller_batches:
    input:
        toxteller_batch_reports
    output:
        report=(
            "results/tox_check/toxteller/{peptide_set}/"
            "clusters_{peptide_set}_rep_seq_toxteller.csv"
        ),
        validated=touch(
            "results/tox_check/toxteller/{peptide_set}/.batches_validated"
        ),
    conda:
        "../envs/tox_check/toxteller.yml"
    shell:
        r"""
        python code/merge_csv_reports.py \
            --inputs {input} \
            --output {output.report}
        """

rule amptox_batch:
    input:
        fasta=tox_check_batch_fasta,
        mapping=tox_check_batch_mapping,
    output:
        report=(
            "results/tox_check/amptox/{peptide_set}/batches/"
            "batch_{batch_id}.csv"
        ),
    params:
        outdir="results/tox_check/amptox/{peptide_set}/batches",
        workdir="results/tox_check/amptox/{peptide_set}/work/batch_{batch_id}",
        raw_report="results/tox_check/amptox/{peptide_set}/work/batch_{batch_id}/raw_amptox.csv",
    threads: 1
    resources:
        mem_mb=get_mem_mb
    conda:
        "../envs/tox_check/amptox.yml"
    shell:
        r"""
        mkdir -p {params.outdir}
        rm -rf {params.workdir}
        mkdir -p {params.workdir}
        input_fasta="$(realpath {input.fasta})"
        raw_report="$(realpath -m {params.raw_report})"
        output_report="$(realpath -m {output.report})"
        python code/tox_check/run_amptox.py \
            --fasta "$input_fasta" \
            --resources resources/amptox \
            --output "$raw_report"
        python code/tox_check/annotate_amptox_report.py \
            --raw-report "$raw_report" \
            --mapping {input.mapping} \
            --output "$output_report"
        """

rule merge_amptox_batches:
    input:
        lambda wildcards: tox_check_batch_reports(wildcards, "amptox")
    output:
        report=(
            "results/tox_check/amptox/{peptide_set}/"
            "clusters_{peptide_set}_rep_seq_amptox.csv"
        ),
        validated=touch(
            "results/tox_check/amptox/{peptide_set}/.batches_validated"
        ),
    conda:
        "../envs/tox_check/amptox.yml"
    shell:
        r"""
        python code/merge_csv_reports.py \
            --inputs {input} \
            --output {output.report}
        """


rule filter_captp_batch:
    input:
        fasta=tox_check_batch_fasta,
    output:
        fasta=(
            f"{CAPTP_BATCH_DIR}/{{peptide_set}}/batch_{{batch_id}}.fasta"
        ),
        stats=(
            f"{CAPTP_BATCH_DIR}/{{peptide_set}}/batch_{{batch_id}}.stats"
        ),
    params:
        max_sequence_length=49,
    conda:
        "../envs/tox_check/toxinpred3_captp.yml"
    shell:
        r"""
        python code/tox_check/filter_fasta_by_length.py \
            --input {input.fasta} \
            --output {output.fasta} \
            --stats {output.stats} \
            --min-length 1 \
            --max-length {params.max_sequence_length}
        """


rule captp_batch:
    input:
        resources_checked=ancient(".snakemake/checks/external_resources_checked"),
        fasta=ancient(
            f"{CAPTP_BATCH_DIR}/{{peptide_set}}/batch_{{batch_id}}.fasta"
        ),
    output:
        report=(
            "results/tox_check/captp/{peptide_set}/batches/"
            "batch_{batch_id}.csv"
        ),
    params:
        outdir="results/tox_check/captp/{peptide_set}/batches",
        tool_dir=CAPTP_PROGRAM_DIR,
        report_name="batch_{batch_id}.csv",
    threads: 1
    resources:
        mem_mb=get_mem_mb
    conda:
        "../envs/tox_check/toxinpred3_captp.yml"
    shell:
        r"""
        mkdir -p {params.outdir}
        input_fasta="$(realpath {input.fasta})"
        output_report="$(realpath -m {output.report})"
        if [ "$(grep -c '^>' "$input_fasta")" -eq 0 ]; then
            printf 'Seq_ID,Sequences,Prediction,Confidence\n' > "$output_report"
            exit 0
        fi
        cd {params.tool_dir}
        rm -f "Results/{params.report_name}"
        python main.py -i "$input_fasta" -o "{params.report_name}"
        mv "Results/{params.report_name}" "$output_report"
        """


rule merge_captp_batches:
    input:
        lambda wildcards: tox_check_batch_reports(wildcards, "captp")
    output:
        report=(
            "results/tox_check/captp/{peptide_set}/"
            "clusters_{peptide_set}_rep_seq_captp.csv"
        ),
        validated=touch(
            "results/tox_check/captp/{peptide_set}/.batches_validated"
        ),
    conda:
        "../envs/tox_check/toxinpred3_captp.yml"
    shell:
        r"""
        python code/merge_csv_reports.py \
            --inputs {input} \
            --output {output.report}
        """


rule build_toxicity_summary:
    input:
        script="code/tox_check/build_toxicity_summary.py",
        fasta=(
            "data/curated_md-lais/mmseqs2/{peptide_set}/"
            "clusters_{peptide_set}_rep_seq.fasta"
        ),
        toxinpred3=ancient(
            "results/tox_check/toxinpred3/{peptide_set}/"
            "clusters_{peptide_set}_rep_seq_toxinpred3.csv"
        ),
        toxteller=ancient(
            "results/tox_check/toxteller/{peptide_set}/"
            "clusters_{peptide_set}_rep_seq_toxteller.csv"
        ),
        captp=ancient(
            "results/tox_check/captp/{peptide_set}/"
            "clusters_{peptide_set}_rep_seq_captp.csv"
        ),
        # amptox=ancient(
        #     "results/tox_check/amptox/{peptide_set}/"
        #     "clusters_{peptide_set}_rep_seq_amptox.csv"
        # ),
    output:
        summary=(
            "results/tox_check/toxicity_summary/{peptide_set}/"
            "clusters_{peptide_set}_toxicity_summary.csv"
        ),
    conda:
        "../envs/tox_check/toxinpred3_captp.yml"
    shell:
        r"""
        python {input.script} \
            --fasta {input.fasta} \
            --toxinpred3 {input.toxinpred3} \
            --toxteller {input.toxteller} \
            --captp {input.captp} \
            # --amptox input.amptox \
            --output-csv {output.summary}
        """


rule filter_non_toxic_fasta:
    input:
        script="code/tox_check/filter_fasta_by_toxicity_summary.py",
        fasta=(
            "data/curated_md-lais/mmseqs2/{peptide_set}/"
            "clusters_{peptide_set}_rep_seq.fasta"
        ),
        summary=ancient(
            "results/tox_check/toxicity_summary/{peptide_set}/"
            "clusters_{peptide_set}_toxicity_summary.csv"
        ),
    output:
        fasta=(
            f"{NON_TOXIC_FASTA_ROOT}/{{peptide_set}}/"
            "clusters_{peptide_set}_rep_seq_non_toxic.fasta"
        ),
        stats=(
            f"{NON_TOXIC_FASTA_ROOT}/{{peptide_set}}/"
            "clusters_{peptide_set}_rep_seq_non_toxic.stats"
        ),
    conda:
        "../envs/tox_check/toxinpred3_captp.yml"
    shell:
        r"""
        python {input.script} \
            --fasta {input.fasta} \
            --summary {input.summary} \
            --output-fasta {output.fasta} \
            --stats {output.stats}
        """


checkpoint split_non_toxic_batches:
    input:
        script="code/split_fasta_batches.py",
        fasta=ancient(
            f"{NON_TOXIC_FASTA_ROOT}/{{peptide_set}}/"
            "clusters_{peptide_set}_rep_seq_non_toxic.fasta"
        ),
    output:
        directory(f"{NON_TOXIC_FASTA_ROOT}/{{peptide_set}}/batches"),
    params:
        n_batches=lambda wildcards: n_batches(wildcards.peptide_set),
        id_prefix=lambda wildcards: f"non_toxic_{wildcards.peptide_set}",
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


def non_toxic_batch_dir(wildcards):
    return checkpoints.split_non_toxic_batches.get(
        peptide_set=wildcards.peptide_set
    ).output[0]


def non_toxic_batch_fasta(wildcards):
    return ancient(
        f"{non_toxic_batch_dir(wildcards)}/batch_{wildcards.batch_id}.fasta"
    )


def non_toxic_batch_mapping(wildcards):
    return ancient(
        f"{non_toxic_batch_dir(wildcards)}/batch_{wildcards.batch_id}.mapping.csv"
    )
