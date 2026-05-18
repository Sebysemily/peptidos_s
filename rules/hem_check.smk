THREADS = int(config.get("max_threads", 8))


rule hemopi2_classification_batch:
    input:
        batch_dir=tox_check_batch_dir,
    output:
        report=(
            "results/hemopi2_classification/{peptide_set}/batches/"
            "batch_{batch_id}.csv"
        ),
    params:
        fasta=lambda wildcards, input: (
            f"{input.batch_dir}/batch_{wildcards.batch_id}.fasta"
        ),
        outdir="results/hemopi2_classification/{peptide_set}/batches",
        workdir="results/hemopi2_classification/{peptide_set}/work/batch_{batch_id}",
        report_name="batch_{batch_id}.csv",
    threads: THREADS
    conda:
        "../envs/hemolisis_check.yml"
    shell:
        r"""
        mkdir -p {params.outdir}
        rm -rf {params.workdir}
        mkdir -p {params.workdir}
        input_fasta="$(realpath {params.fasta})"
        output_report="$(realpath -m {output.report})"
        workdir="$(realpath -m {params.workdir})"

        if command -v hemopi2_classification >/dev/null 2>&1; then
            hemopi2_cmd="hemopi2_classification"
        elif command -v hemopi2_classification.py >/dev/null 2>&1; then
            hemopi2_cmd="hemopi2_classification.py"
        else
            echo "HemoPI2 classification command not found." >&2
            echo "Expected hemopi2_classification or hemopi2_classification.py." >&2
            exit 1
        fi

        "$hemopi2_cmd" \
            -i "$input_fasta" \
            -o "{params.report_name}" \
            -j 1 \
            -m 2 \
            -d 2 \
            -wd "$workdir"

        for candidate in \
            "$workdir/{params.report_name}" \
            "$workdir/final_output.csv" \
            "$workdir/final_output" \
            "$workdir/outfile.csv"; do
            if [ -f "$candidate" ]; then
                mv "$candidate" "$output_report"
                break
            fi
        done

        if [ ! -f "$output_report" ]; then
            echo "HemoPI2 classification did not create an expected CSV output." >&2
            exit 1
        fi
        """


rule merge_hemopi2_classification_batches:
    input:
        lambda wildcards: tox_check_batch_reports(
            wildcards, "hemopi2_classification"
        )
    output:
        report=(
            "results/hemopi2_classification/{peptide_set}/"
            "clusters_{peptide_set}_rep_seq_hemopi2_classification.csv"
        ),
    conda:
        "../envs/tox_check.yml"
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
            "results/hemopi2_regression/{peptide_set}/batches/"
            "batch_{batch_id}.csv"
        ),
    params:
        fasta=lambda wildcards, input: (
            f"{input.batch_dir}/batch_{wildcards.batch_id}.fasta"
        ),
        outdir="results/hemopi2_regression/{peptide_set}/batches",
        workdir="results/hemopi2_regression/{peptide_set}/work/batch_{batch_id}",
        report_name="batch_{batch_id}.csv",
    threads: THREADS
    conda:
        "../envs/hemolisis_check.yml"
    shell:
        r"""
        mkdir -p {params.outdir}
        rm -rf {params.workdir}
        mkdir -p {params.workdir}
        input_fasta="$(realpath {params.fasta})"
        output_report="$(realpath -m {output.report})"
        workdir="$(realpath -m {params.workdir})"

        if command -v hemopi2_regression >/dev/null 2>&1; then
            hemopi2_cmd="hemopi2_regression"
        elif command -v hemopi2_regression.py >/dev/null 2>&1; then
            hemopi2_cmd="hemopi2_regression.py"
        else
            echo "HemoPI2 regression command not found." >&2
            echo "Expected hemopi2_regression or hemopi2_regression.py." >&2
            exit 1
        fi

        "$hemopi2_cmd" \
            -i "$input_fasta" \
            -o "{params.report_name}" \
            -j 1 \
            -d 2 \
            -wd "$workdir"

        for candidate in \
            "$workdir/{params.report_name}" \
            "$workdir/final_output.csv" \
            "$workdir/final_output" \
            "$workdir/outfile.csv"; do
            if [ -f "$candidate" ]; then
                mv "$candidate" "$output_report"
                break
            fi
        done

        if [ ! -f "$output_report" ]; then
            echo "HemoPI2 regression did not create an expected CSV output." >&2
            exit 1
        fi
        """


rule merge_hemopi2_regression_batches:
    input:
        lambda wildcards: tox_check_batch_reports(
            wildcards, "hemopi2_regression"
        )
    output:
        report=(
            "results/hemopi2_regression/{peptide_set}/"
            "clusters_{peptide_set}_rep_seq_hemopi2_regression.csv"
        ),
    conda:
        "../envs/tox_check.yml"
    shell:
        r"""
        python code/merge_csv_reports.py \
            --inputs {input} \
            --output {output.report}
        """
