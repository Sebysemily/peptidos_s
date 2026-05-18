import re

THREADS = int(config.get("max_threads", 8))

PEPTIDE_SET = "|".join(
    re.escape(name) for name in config["curated_fastas"]
)


TOXTELLER_CONFIG = config.get("toxteller", {})
TOXTELLER_PROGRAM_DIR = TOXTELLER_CONFIG["program_dir"]



wildcard_constraints:
    peptide_set=PEPTIDE_SET,
    batch_id=r"\d+"


rule split_toxinpred3_batches:
    input:
        rep_seq=(
            "data/curated_md-lais/mmseqs2/{peptide_set}/"
            "clusters_{peptide_set}_rep_seq.fasta"
        ),
    output:
        complete="results/batches/toxinpred3/{peptide_set}/.complete",
    params:
        outdir="results/batches/toxinpred3/{peptide_set}",
        n_batches=lambda wildcards: n_batches(wildcards.peptide_set),
    conda:
        "../envs/tox_check.yml"
    shell:
        r"""
        python scripts/split_fasta_batches.py \
            --input {input.rep_seq} \
            --outdir {params.outdir} \
            --n-batches {params.n_batches}
        """


rule toxinpred3_batch:
    input:
        batches_done="results/batches/toxinpred3/{peptide_set}/.complete",
    output:
        report=(
            "results/toxinpred3/{peptide_set}/batches/"
            "batch_{batch_id}.csv"
        ),
    params:
        fasta="results/batches/toxinpred3/{peptide_set}/batch_{batch_id}.fasta",
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
        lambda wildcards: expand(
            "results/toxinpred3/{peptide_set}/batches/batch_{batch_id}.csv",
            peptide_set=wildcards.peptide_set,
            batch_id=batch_ids_for(wildcards.peptide_set),
        )
    output:
        report=(
            "results/toxinpred3/{peptide_set}/"
            "clusters_{peptide_set}_rep_seq_toxinpred3.csv"
        ),
    conda:
        "../envs/tox_check.yml"
    shell:
        r"""
        python scripts/merge_csv_reports.py \
            --inputs {input} \
            --output {output.report}
        """


checkpoint split_toxteller_batches:
    input:
        rep_seq=(
            "data/curated_md-lais/mmseqs2/{peptide_set}/"
            "clusters_{peptide_set}_rep_seq.fasta"
        ),
    output:
        directory("results/batches/toxteller/{peptide_set}"),
    params:
        max_sequences_per_batch=10000,
    conda:
        "../envs/tox_check.yml"
    shell:
        r"""
        python scripts/split_fasta_batches.py \
            --input {input.rep_seq} \
            --outdir {output} \
            --max-records-per-batch {params.max_sequences_per_batch}
        """


rule toxteller_batch:
    input:
        fasta="results/batches/toxteller/{peptide_set}/batch_{batch_id}.fasta",
    output:
        report=(
            "results/toxteller/{peptide_set}/batches/"
            "batch_{batch_id}.csv"
        ),
    params:
        outdir="results/toxteller/{peptide_set}/batches",
        program_dir=TOXTELLER_PROGRAM_DIR,
        repo_url="https://github.com/comics-asiis/ToxicPeptidePrediction.git",
    threads: THREADS
    conda:
        "../envs/tox_check.yml"
    shell:
        r"""
        mkdir -p {params.outdir}
        input_fasta="$(realpath {input.fasta})"
        output_report="$(realpath -m {output.report})"
        if [ ! -f "{params.program_dir}/program_resource/toxteller.py" ]; then
            if [ -d "{params.program_dir}" ] && [ -n "$(find "{params.program_dir}" -mindepth 1 -maxdepth 1 -print -quit)" ]; then
                echo "ToxTeller checkout is incomplete: {params.program_dir}" >&2
                echo "Remove it or set toxteller.program_dir to a valid ToxTeller checkout." >&2
                exit 1
            fi
            mkdir -p "$(dirname "{params.program_dir}")"
            git clone "{params.repo_url}" "{params.program_dir}"
        fi
        cd {params.program_dir}/program_resource
        python toxteller.py "$input_fasta"
        mv "${{input_fasta}}.csv" "$output_report"
        """


def toxteller_batch_reports(wildcards):
    batch_dir = checkpoints.split_toxteller_batches.get(
        peptide_set=wildcards.peptide_set
    ).output[0]
    batch_ids = glob_wildcards(
        f"{batch_dir}/batch_{{batch_id}}.fasta"
    ).batch_id
    batch_ids = sorted(batch_ids, key=int)
    return expand(
        "results/toxteller/{peptide_set}/batches/batch_{batch_id}.csv",
        peptide_set=wildcards.peptide_set,
        batch_id=batch_ids,
    )


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
        python scripts/merge_csv_reports.py \
            --inputs {input} \
            --output {output.report}
        """
