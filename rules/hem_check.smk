import re

THREADS = int(config.get("max_threads", 8))
HEMOPI2_THREADS = 1
MACREL_THREADS = THREADS
HEPAD_THREADS = 1
HEMO_DL_THREADS = THREADS
HEPAD_PROGRAM_DIR = "resources/HEPAD"
VALID_HEPAD_DATASETS = {"Hmp1", "Hmp2", "Hmp3", "Hmpm"}
HEPAD_DATASETS = config.get(
    "hepad_datasets", ["Hmp1", "Hmp2", "Hmp3", "Hmpm"]
)
if not HEPAD_DATASETS:
    raise ValueError("hepad_datasets must contain at least one dataset")
for dataset in HEPAD_DATASETS:
    if dataset not in VALID_HEPAD_DATASETS:
        raise ValueError(
            f"Unknown HEPAD dataset {dataset!r}; "
            f"expected one of {sorted(VALID_HEPAD_DATASETS)}"
        )
HEPAD_DATASET_PATTERN = "|".join(re.escape(dataset) for dataset in HEPAD_DATASETS)


wildcard_constraints:
    hepad_dataset=HEPAD_DATASET_PATTERN,
    batch_id=r"\d+"


def hem_check_batch_reports(wildcards, tool):
    batch_dir = non_toxic_batch_dir(wildcards)
    batch_ids = glob_wildcards(
        f"{batch_dir}/batch_{{batch_id}}.fasta"
    ).batch_id
    batch_ids = sorted(batch_ids, key=int)
    reports = expand(
        "results/hemo_check/{tool}/{peptide_set}/batches/batch_{batch_id}.csv",
        tool=tool,
        peptide_set=wildcards.peptide_set,
        batch_id=batch_ids,
    )
    return [ancient(report) for report in reports]


def hepad_batch_reports(wildcards):
    batch_dir = non_toxic_batch_dir(wildcards)
    batch_ids = glob_wildcards(
        f"{batch_dir}/batch_{{batch_id}}.fasta"
    ).batch_id
    batch_ids = sorted(batch_ids, key=int)
    reports = expand(
        (
            "results/hemo_check/hepad/{hepad_dataset}/{peptide_set}/"
            "batches/batch_{batch_id}.csv"
        ),
        hepad_dataset=wildcards.hepad_dataset,
        peptide_set=wildcards.peptide_set,
        batch_id=batch_ids,
    )
    return [ancient(report) for report in reports]


rule hemopi2_classification_batch:
    input:
        fasta=non_toxic_batch_fasta,
        mapping=non_toxic_batch_mapping,
    output:
        report=(
            "results/hemo_check/hemopi2_classification/{peptide_set}/batches/"
            "batch_{batch_id}.csv"
        ),
    params:
        workdir="results/hemo_check/hemopi2_classification/{peptide_set}/work/batch_{batch_id}",
        raw_report="results/hemo_check/hemopi2_classification/{peptide_set}/work/batch_{batch_id}/raw_hemopi2_classification.csv",
        raw_report_name="raw_hemopi2_classification.csv",
    threads: HEMOPI2_THREADS
    resources:
        mem_mb=get_mem_mb
    conda:
        "../envs/hem_check/hemopi2.yml"
    shell:
        r"""
        rm -rf {params.workdir}
        mkdir -p {params.workdir}
        output_report="$(realpath -m {output.report})"
        workdir="$(realpath -m {params.workdir})"
        mkdir -p "$(dirname "$output_report")"

        input_fasta="$(realpath {input.fasta})"

        hemopi2_classification \
            -i "$input_fasta" \
            -o {params.raw_report_name} \
            -j 1 \
            -m 2 \
            -d 2 \
            -wd "$workdir"

        python code/hem_check/annotate_indexed_report.py \
            --raw-report {params.raw_report} \
            --mapping {input.mapping} \
            --output "$output_report" \
            --prefix hemopi2_classification
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
        fasta=non_toxic_batch_fasta,
        mapping=non_toxic_batch_mapping,
    output:
        report=(
            "results/hemo_check/hemopi2_regression/{peptide_set}/batches/"
            "batch_{batch_id}.csv"
        ),
    params:
        workdir="results/hemo_check/hemopi2_regression/{peptide_set}/work/batch_{batch_id}",
        raw_report="results/hemo_check/hemopi2_regression/{peptide_set}/work/batch_{batch_id}/raw_hemopi2_regression.csv",
        raw_report_name="raw_hemopi2_regression.csv",
    threads: HEMOPI2_THREADS
    resources:
        mem_mb=get_mem_mb
    conda:
        "../envs/hem_check/hemopi2.yml"
    shell:
        r"""
        rm -rf {params.workdir}
        mkdir -p {params.workdir}
        output_report="$(realpath -m {output.report})"
        workdir="$(realpath -m {params.workdir})"
        mkdir -p "$(dirname "$output_report")"

        input_fasta="$(realpath {input.fasta})"

        hemopi2_regression \
            -i "$input_fasta" \
            -o {params.raw_report_name} \
            -j 1 \
            -d 2 \
            -wd "$workdir"

        python code/hem_check/annotate_indexed_report.py \
            --raw-report {params.raw_report} \
            --mapping {input.mapping} \
            --output "$output_report" \
            --prefix hemopi2_regression
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


rule macrel_batch:
    input:
        fasta=non_toxic_batch_fasta,
        mapping=non_toxic_batch_mapping,
    output:
        report=(
            "results/hemo_check/macrel/{peptide_set}/batches/"
            "batch_{batch_id}.csv"
        ),
    params:
        outdir="results/hemo_check/macrel/{peptide_set}/batches",
        workdir="results/hemo_check/macrel/{peptide_set}/work/batch_{batch_id}",
        tag="batch_{batch_id}",
    threads: MACREL_THREADS
    resources:
        mem_mb=get_mem_mb
    conda:
        "../envs/hem_check/macrel.yml"
    shell:
        r"""
        mkdir -p {params.outdir}
        rm -rf {params.workdir}
        mkdir -p {params.workdir}

        input_fasta="$(realpath {input.fasta})"
        workdir="$(realpath -m {params.workdir})"

        macrel peptides \
            --fasta "$input_fasta" \
            --output "$workdir" \
            --threads {threads} \
            --tag {params.tag} \
            -- query-mode=mmseqs \
            -- local \
            --keep-negatives \
            --force

        python code/hem_check/annotate_macrel_report.py \
            --input-dir "$workdir" \
            --mapping {input.mapping} \
            --output {output.report}
        """


rule merge_macrel_batches:
    input:
        lambda wildcards: hem_check_batch_reports(wildcards, "macrel")
    output:
        report=(
            "results/hemo_check/macrel/{peptide_set}/"
            "clusters_{peptide_set}_rep_seq_macrel.csv"
        ),
    conda:
        "../envs/tox_check/toxinpred3_captp.yml"
    shell:
        r"""
        python code/merge_csv_reports.py \
            --inputs {input} \
            --output {output.report}
        """


rule hepad_batch:
    input:
        fasta=non_toxic_batch_fasta,
        mapping=non_toxic_batch_mapping,
    output:
        report=(
            "results/hemo_check/hepad/{hepad_dataset}/{peptide_set}/batches/"
            "batch_{batch_id}.csv"
        ),
    params:
        hepad_root=HEPAD_PROGRAM_DIR,
        workdir=(
            "results/hemo_check/hepad/{hepad_dataset}/{peptide_set}/"
            "work/batch_{batch_id}"
        ),
        binary_report=(
            "results/hemo_check/hepad/{hepad_dataset}/{peptide_set}/"
            "work/batch_{batch_id}/output/binary_vector.csv"
        ),
        probability_report=(
            "results/hemo_check/hepad/{hepad_dataset}/{peptide_set}/"
            "work/batch_{batch_id}/output/probability_vector.csv"
        ),
        prefix="hepad_{hepad_dataset}",
    threads: HEPAD_THREADS
    resources:
        mem_mb=get_hepad_mem_mb
    conda:
        "../envs/hem_check/hepad.yml"
    shell:
        r"""
        rm -rf {params.workdir}
        mkdir -p {params.workdir}
        output_report="$(realpath -m {output.report})"
        workdir="$(realpath -m {params.workdir})"
        mkdir -p "$(dirname "$output_report")"

        input_fasta="$(realpath {input.fasta})"

        python -c "import modlamp" 2>/dev/null || python -m pip install --no-deps modlamp==4.3.0

        python code/hem_check/run_hepad.py \
            --input "$input_fasta" \
            --dataset {wildcards.hepad_dataset} \
            --hepad-root {params.hepad_root} \
            --workdir "$workdir"

        python code/hem_check/annotate_hepad_report.py \
            --binary-report {params.binary_report} \
            --probability-report {params.probability_report} \
            --mapping {input.mapping} \
            --output "$output_report" \
            --prefix {params.prefix}
        """


rule merge_hepad_batches:
    input:
        hepad_batch_reports
    output:
        report=(
            "results/hemo_check/hepad/{hepad_dataset}/{peptide_set}/"
            "clusters_{peptide_set}_rep_seq_hepad_{hepad_dataset}.csv"
        ),
    conda:
        "../envs/tox_check/toxinpred3_captp.yml"
    shell:
        r"""
        python code/merge_csv_reports.py \
            --inputs {input} \
            --output {output.report}
        """


rule hemo_DL_batch:
    input:
        fasta=non_toxic_batch_fasta,
        mapping=non_toxic_batch_mapping,
    output:
        report=(
            "results/hemo_check/hemo_DL/{peptide_set}/batches/"
            "batch_{batch_id}.csv"
        ),
    params:
        workdir="results/hemo_check/hemo_DL/{peptide_set}/work/batch_{batch_id}",
        script_dir="resources/HemoDL/source",
        raw_report="results/hemo_check/hemo_DL/{peptide_set}/work/batch_{batch_id}/predict_results.csv",
        raw_report_header="results/hemo_check/hemo_DL/{peptide_set}/work/batch_{batch_id}/predict_results_header.csv",
    threads: HEMO_DL_THREADS
    resources:
        mem_mb=get_mem_mb
    conda:
        "../envs/hem_check/hemo_DL.yml"
    shell:
        r"""
        rm -rf {params.workdir}
        mkdir -p {params.workdir}
        output_report="$(realpath -m {output.report})"
        workdir="$(realpath -m {params.workdir})"
        mkdir -p "$(dirname "$output_report")"

        input_fasta="$(realpath {input.fasta})"
        script_dir_abs="$(realpath {params.script_dir})"

        # Copy models and source files to avoid parallel write conflicts on predict_results.csv
        cp -r "$script_dir_abs/models" "$workdir/"
        cp "$script_dir_abs/features.py" "$workdir/"
        cp "$script_dir_abs/predict.py" "$workdir/"

        cd "$workdir"
        python predict.py -p "$input_fasta"
        cd -

        # Add header to the raw report since predict.py output doesn't contain one
        echo "id,probability" > {params.raw_report_header}
        cat {params.raw_report} >> {params.raw_report_header}

        python code/hem_check/annotate_indexed_report.py \
            --raw-report {params.raw_report_header} \
            --mapping {input.mapping} \
            --output "$output_report" \
            --prefix hemo_DL
        """


rule build_hemolytic_summary:
    input:
        script="code/hem_check/build_hemolytic_summary.py",
        fasta=(
            f"{NON_TOXIC_FASTA_ROOT}/{{peptide_set}}/"
            "clusters_{peptide_set}_rep_seq_non_toxic.fasta"
        ),
        macrel=ancient(
            "results/hemo_check/macrel/{peptide_set}/"
            "clusters_{peptide_set}_rep_seq_macrel.csv"
        ),
        hemopi2_classification=ancient(
            "results/hemo_check/hemopi2_classification/{peptide_set}/"
            "clusters_{peptide_set}_rep_seq_hemopi2_classification.csv"
        ),
        hemopi2_regression=ancient(
            "results/hemo_check/hemopi2_regression/{peptide_set}/"
            "clusters_{peptide_set}_rep_seq_hemopi2_regression.csv"
        ),
        hepad_hmp1=ancient(
            "results/hemo_check/hepad/Hmp1/{peptide_set}/"
            "clusters_{peptide_set}_rep_seq_hepad_Hmp1.csv"
        ),
        hepad_hmpm=ancient(
            "results/hemo_check/hepad/Hmpm/{peptide_set}/"
            "clusters_{peptide_set}_rep_seq_hepad_Hmpm.csv"
        ),
    output:
        summary=(
            "results/hemo_check/hemolytic_summary/{peptide_set}/"
            "clusters_{peptide_set}_hemolytic_summary.csv"
        ),
    conda:
        "../envs/tox_check/toxinpred3_captp.yml"
    shell:
        r"""
        python {input.script} \
            --fasta {input.fasta} \
            --macrel {input.macrel} \
            --hemopi2-classification {input.hemopi2_classification} \
            --hemopi2-regression {input.hemopi2_regression} \
            --hepad-hmp1 {input.hepad_hmp1} \
            --hepad-hmpm {input.hepad_hmpm} \
            --output-csv {output.summary}
        """
rule merge_hemo_DL_batches:
    input:
        lambda wildcards: hem_check_batch_reports(wildcards, "hemo_DL")
    output:
        report=(
            "results/hemo_check/hemo_DL/{peptide_set}/"
            "clusters_{peptide_set}_rep_seq_hemo_DL.csv"
        ),
    conda:
        "../envs/tox_check/toxinpred3_captp.yml"
    shell:
        r"""
        python code/merge_csv_reports.py \
            --inputs {input} \
            --output {output.report}
        """


