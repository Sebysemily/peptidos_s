THREADS = int(config.get("max_threads", 8))


def hem_check_batch_reports(wildcards, tool):
    batch_dir = tox_check_batch_dir(wildcards)
    batch_ids = glob_wildcards(
        f"{batch_dir}/batch_{{batch_id}}.fasta"
    ).batch_id
    batch_ids = sorted(batch_ids, key=int)
    return expand(
        "results/hemo_check/{tool}/{peptide_set}/batches/batch_{batch_id}.csv",
        tool=tool,
        peptide_set=wildcards.peptide_set,
        batch_id=batch_ids,
    )


rule hemopi2_classification_batch:
    input:
        batch_dir=tox_check_batch_dir,
    output:
        report=(
            "results/hemo_check/hemopi2_classification/{peptide_set}/batches/"
            "batch_{batch_id}.csv"
        ),
    params:
        fasta=lambda wildcards, input: (
            f"{input.batch_dir}/batch_{wildcards.batch_id}.fasta"
        ),
        workdir="results/hemo_check/hemopi2_classification/{peptide_set}/work/batch_{batch_id}",
        input_fasta="results/hemo_check/hemopi2_classification/{peptide_set}/work/batch_{batch_id}/input.fasta",
        report_name="batch_{batch_id}.csv",
    threads: THREADS
    conda:
        "../envs/hemolisis_check.yml"
    shell:
        r"""
        rm -rf {params.workdir}
        mkdir -p {params.workdir}
        output_report="$(realpath -m {output.report})"
        workdir="$(realpath -m {params.workdir})"
        mkdir -p "$(dirname "$output_report")"

        python code/hem_check/sanitize_fasta_for_hemopi2.py \
            --input {params.fasta} \
            --output {params.input_fasta}
        input_fasta="$(realpath {params.input_fasta})"

        hemopi2_classification \
            -i "$input_fasta" \
            -o {params.report_name} \
            -j 1 \
            -m 2 \
            -d 2 \
            -wd "$workdir"

        mv "$workdir/{params.report_name}" "$output_report"
        """


rule merge_hemopi2_classification_batches:
    input:
        lambda wildcards: hem_check_batch_reports(
            wildcards, "hemopi2_classification"
        )
    output:
        report=(
            "results/hemo_check/hemopi2_classification/{peptide_set}/"
            "clusters_{peptide_set}_rep_seq_hemopi2_classification.csv"
        ),
    conda:
        "../envs/tox_check/toxinpred3_captp.yml"
    shell:
        r"""
        python code/merge_csv_reports.py \
            --inputs {input} \
            --output {output.report}
        """


rule hemopi2_regression_batch:
    input:
        batch_dir=tox_check_batch_dir,
    output:
        report=(
            "results/hemo_check/hemopi2_regression/{peptide_set}/batches/"
            "batch_{batch_id}.csv"
        ),
    params:
        fasta=lambda wildcards, input: (
            f"{input.batch_dir}/batch_{wildcards.batch_id}.fasta"
        ),
        workdir="results/hemo_check/hemopi2_regression/{peptide_set}/work/batch_{batch_id}",
        input_fasta="results/hemo_check/hemopi2_regression/{peptide_set}/work/batch_{batch_id}/input.fasta",
        report_name="batch_{batch_id}.csv",
    threads: THREADS
    conda:
        "../envs/hemolisis_check.yml"
    shell:
        r"""
        rm -rf {params.workdir}
        mkdir -p {params.workdir}
        output_report="$(realpath -m {output.report})"
        workdir="$(realpath -m {params.workdir})"
        mkdir -p "$(dirname "$output_report")"

        python code/hem_check/sanitize_fasta_for_hemopi2.py \
            --input {params.fasta} \
            --output {params.input_fasta}
        input_fasta="$(realpath {params.input_fasta})"

        hemopi2_regression \
            -i "$input_fasta" \
            -o {params.report_name} \
            -j 1 \
            -d 2 \
            -wd "$workdir"

        mv "$workdir/{params.report_name}" "$output_report"
        """


rule merge_hemopi2_regression_batches:
    input:
        lambda wildcards: hem_check_batch_reports(
            wildcards, "hemopi2_regression"
        )
    output:
        report=(
            "results/hemo_check/hemopi2_regression/{peptide_set}/"
            "clusters_{peptide_set}_rep_seq_hemopi2_regression.csv"
        ),
    conda:
        "../envs/tox_check/toxinpred3_captp.yml"
    shell:
        r"""
        python code/merge_csv_reports.py \
            --inputs {input} \
            --output {output.report}
        """
