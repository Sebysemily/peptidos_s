import re

THREADS = int(config.get("max_threads", 8))

PEPTIDE_SET = "|".join(
    re.escape(name) for name in config["curated_fastas"]
)


TOXTELLER_CONFIG = config.get("toxteller", {})
TOXTELLER_PROGRAM_DIR = TOXTELLER_CONFIG["program_dir"]

CAPTP_CONFIG = config.get("captp", {})
CAPTP_PROGRAM_DIR = CAPTP_CONFIG["program_dir"]


wildcard_constraints:
    peptide_set=PEPTIDE_SET,
    batch_id=r"\d+"


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


checkpoint split_tox_check_batches:
    input:
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
        python code/tox_check/split_fasta_batches.py \
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
        python code/tox_check/split_fasta_batches.py \
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
    threads: THREADS
    conda:
        "../envs/tox_check.yml"
    shell:
        r"""
        mkdir -p {params.outdir}
        rm -rf {params.workdir}
        mkdir -p {params.workdir}
        input_fasta="$(realpath {params.fasta})"
        output_report="$(realpath -m {output.report})"
        cd {params.workdir}
        toxinpred3 -i "$input_fasta" -o "$output_report" -m 2 -d 2
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
        program_dir=TOXTELLER_PROGRAM_DIR,
    threads: THREADS
    conda:
        "../envs/tox_check.yml"
    shell:
        r"""
        mkdir -p {params.outdir}
        input_fasta="$(realpath {params.fasta})"
        output_report="$(realpath -m {output.report})"
        if [ ! -f "{params.program_dir}/program_resource/toxteller.py" ]; then
            echo "ToxTeller checkout is missing or incomplete: {params.program_dir}" >&2
            echo "Initialize resources with: git submodule update --init --recursive" >&2
            echo "Expected file: {params.program_dir}/program_resource/toxteller.py" >&2
            exit 1
        fi
        cd {params.program_dir}/program_resource
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
        program_dir=CAPTP_PROGRAM_DIR,
        report_name="batch_{batch_id}.csv",
    threads: THREADS
    conda:
        "../envs/tox_check.yml"
    shell:
        r"""
        mkdir -p {params.outdir}
        input_fasta="$(realpath {input.fasta})"
        output_report="$(realpath -m {output.report})"
        if [ ! -f "{params.program_dir}/main.py" ]; then
            echo "CAPTP checkout is missing or incomplete: {params.program_dir}" >&2
            echo "Initialize resources with: git submodule update --init --recursive" >&2
            echo "Expected file: {params.program_dir}/main.py" >&2
            exit 1
        fi
        if [ "$(grep -c '^>' "$input_fasta")" -eq 0 ]; then
            printf 'Seq_ID,Sequences,Prediction,Confidence\n' > "$output_report"
            exit 0
        fi
        cd {params.program_dir}
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
