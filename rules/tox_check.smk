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
        script="code/tox_check/split_fasta_batches.py",
        rep_seq=(
            "data/curated_md-lais/mmseqs2/{peptide_set}/"
            "clusters_{peptide_set}_rep_seq.fasta"
        ),
    output:
        directory("results/batches/tox_check/{peptide_set}"),
    params:
        n_batches=lambda wildcards: n_batches(wildcards.peptide_set),
    conda:
        "../envs/tox_check.yml"
    shell:
        r"""
        python {input.script} \
            --input {input.rep_seq} \
            --outdir {output} \
            --n-batches {params.n_batches}
        """


def tox_check_batch_dir(wildcards):
    return checkpoints.split_tox_check_batches.get(
        peptide_set=wildcards.peptide_set
    ).output[0]


def tox_check_batch_reports(wildcards, tool):
    batch_dir = tox_check_batch_dir(wildcards)
    batch_ids = glob_wildcards(
        f"{batch_dir}/batch_{{batch_id}}.fasta"
    ).batch_id
    batch_ids = sorted(batch_ids, key=int)
    return expand(
        "results/{tool}/{peptide_set}/batches/batch_{batch_id}.csv",
        tool=tool,
        peptide_set=wildcards.peptide_set,
        batch_id=batch_ids,
    )


checkpoint split_toxteller_batches:
    input:
        script="code/tox_check/split_fasta_batches.py",
        rep_seq=(
            "data/curated_md-lais/mmseqs2/{peptide_set}/"
            "clusters_{peptide_set}_rep_seq.fasta"
        ),
    output:
        directory("results/batches/toxteller/{peptide_set}"),
    params:
        max_sequences_per_batch=9500,
    conda:
        "../envs/tox_check.yml"
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


def toxteller_batch_reports(wildcards):
    batch_dir = toxteller_batch_dir(wildcards)
    batch_ids = glob_wildcards(
        f"{batch_dir}/batch_{{batch_id}}.fasta"
    ).batch_id
    batch_ids = sorted(batch_ids, key=int)
    return expand(
        "results/toxteller/{peptide_set}/batches/batch_{batch_id}.csv",
        peptide_set=wildcards.peptide_set,
        batch_id=batch_ids,
    )


rule toxinpred3_batch:
    input:
        batch_dir=tox_check_batch_dir,
    output:
        report=(
            "results/toxinpred3/{peptide_set}/batches/"
            "batch_{batch_id}.csv"
        ),
    params:
        fasta=lambda wildcards, input: (
            f"{input.batch_dir}/batch_{wildcards.batch_id}.fasta"
        ),
        outdir="results/toxinpred3/{peptide_set}/batches",
        workdir="results/toxinpred3/{peptide_set}/work/batch_{batch_id}",
        indexed_fasta="results/toxinpred3/{peptide_set}/work/batch_{batch_id}/input_indexed.fasta",
        mapping="results/toxinpred3/{peptide_set}/work/batch_{batch_id}/input_mapping.csv",
        raw_report="results/toxinpred3/{peptide_set}/work/batch_{batch_id}/raw_toxinpred3.csv",
    threads: THREADS
    conda:
        "../envs/tox_check.yml"
    shell:
        r"""
        mkdir -p {params.outdir}
        rm -rf {params.workdir}
        mkdir -p {params.workdir}
        python code/tox_check/prepare_toxinpred3_input.py \
            --input {params.fasta} \
            --output-fasta {params.indexed_fasta} \
            --mapping {params.mapping} \
            --prefix "{wildcards.peptide_set}_{wildcards.batch_id}"
        input_fasta="$(realpath {params.indexed_fasta})"
        raw_report="$(realpath -m {params.raw_report})"
        output_report="$(realpath -m {output.report})"
        cd {params.workdir}
        toxinpred3 -i "$input_fasta" -o "$raw_report" -m 2 -d 2
        cd -
        python code/tox_check/annotate_toxinpred3_report.py \
            --raw-report "$raw_report" \
            --mapping {params.mapping} \
            --output "$output_report"
        """


rule merge_toxinpred3_batches:
    input:
        lambda wildcards: tox_check_batch_reports(wildcards, "toxinpred3")
    output:
        report=(
            "results/toxinpred3/{peptide_set}/"
            "clusters_{peptide_set}_rep_seq_toxinpred3.csv"
        ),
    conda:
        "../envs/tox_check.yml"
    shell:
        r"""
        python code/merge_csv_reports.py \
            --inputs {input} \
            --output {output.report}
        """

rule toxteller_batch:
    input:
        resources_checked="results/setup/.external_resources_checked",
        batch_dir=toxteller_batch_dir,
    output:
        report=(
            "results/toxteller/{peptide_set}/batches/"
            "batch_{batch_id}.csv"
        ),
    params:
        fasta=lambda wildcards, input: (
            f"{input.batch_dir}/batch_{wildcards.batch_id}.fasta"
        ),
        outdir="results/toxteller/{peptide_set}/batches",
        tool_dir=TOXTELLER_PROGRAM_DIR,
    threads: THREADS
    conda:
        "../envs/tox_check.yml"
    shell:
        r"""
        mkdir -p {params.outdir}
        input_fasta="$(realpath {params.fasta})"
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
            "results/toxteller/{peptide_set}/"
            "clusters_{peptide_set}_rep_seq_toxteller.csv"
        ),
    conda:
        "../envs/tox_check.yml"
    shell:
        r"""
        python code/merge_csv_reports.py \
            --inputs {input} \
            --output {output.report}
        """


rule filter_captp_batch:
    input:
        batch_dir=tox_check_batch_dir,
    output:
        fasta=(
            "results/batches/captp/{peptide_set}/"
            "batch_{batch_id}.fasta"
        ),
        stats=(
            "results/batches/captp/{peptide_set}/"
            "batch_{batch_id}.stats"
        ),
    params:
        fasta=lambda wildcards, input: (
            f"{input.batch_dir}/batch_{wildcards.batch_id}.fasta"
        ),
        max_sequence_length=49,
    conda:
        "../envs/tox_check.yml"
    shell:
        r"""
        python code/tox_check/filter_fasta_by_length.py \
            --input {params.fasta} \
            --output {output.fasta} \
            --stats {output.stats} \
            --min-length 1 \
            --max-length {params.max_sequence_length}
        """


rule captp_batch:
    input:
        resources_checked="results/setup/.external_resources_checked",
        fasta="results/batches/captp/{peptide_set}/batch_{batch_id}.fasta",
    output:
        report=(
            "results/captp/{peptide_set}/batches/"
            "batch_{batch_id}.csv"
        ),
    params:
        outdir="results/captp/{peptide_set}/batches",
        tool_dir=CAPTP_PROGRAM_DIR,
        report_name="batch_{batch_id}.csv",
    threads: THREADS
    conda:
        "../envs/tox_check.yml"
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
            "results/captp/{peptide_set}/"
            "clusters_{peptide_set}_rep_seq_captp.csv"
        ),
    conda:
        "../envs/tox_check.yml"
    shell:
        r"""
        python code/merge_csv_reports.py \
            --inputs {input} \
            --output {output.report}
        """


rule build_toxicity_summary:
    input:
        fasta=(
            "data/curated_md-lais/mmseqs2/{peptide_set}/"
            "clusters_{peptide_set}_rep_seq.fasta"
        ),
        toxinpred3=(
            "results/toxinpred3/{peptide_set}/"
            "clusters_{peptide_set}_rep_seq_toxinpred3.csv"
        ),
        toxteller=(
            "results/toxteller/{peptide_set}/"
            "clusters_{peptide_set}_rep_seq_toxteller.csv"
        ),
        captp=(
            "results/captp/{peptide_set}/"
            "clusters_{peptide_set}_rep_seq_captp.csv"
        ),
    output:
        summary=(
            "results/toxicity_summary/{peptide_set}/"
            "clusters_{peptide_set}_toxicity_summary.csv"
        ),
    conda:
        "../envs/tox_check.yml"
    shell:
        r"""
        python code/tox_check/build_toxicity_summary.py \
            --fasta {input.fasta} \
            --toxinpred3 {input.toxinpred3} \
            --toxteller {input.toxteller} \
            --captp {input.captp} \
            --output-csv {output.summary}
        """
