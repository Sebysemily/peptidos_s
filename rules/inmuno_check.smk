import re

wildcard_constraints:
    peptide_set="|".join(re.escape(name) for name in config["curated_fastas"]),
    batch_id=r"\d+"

def non_hemo_batch_dir(wildcards):
    return checkpoints.split_non_hemo_batches.get(
        peptide_set=wildcards.peptide_set
    ).output[0]

def non_hemo_batch_fasta(wildcards):
    return ancient(f"{non_hemo_batch_dir(wildcards)}/batch_{wildcards.batch_id}.fasta")

def non_hemo_batch_mapping(wildcards):
    return ancient(f"{non_hemo_batch_dir(wildcards)}/batch_{wildcards.batch_id}.mapping.csv")

def inmuno_check_batch_reports(wildcards, tool):
    batch_dir = non_hemo_batch_dir(wildcards)
    batch_ids = glob_wildcards(
        f"{batch_dir}/batch_{{batch_id}}.fasta"
    ).batch_id
    batch_ids = sorted(batch_ids, key=int)
    reports = expand(
        "results/inmuno_check/{tool}/{peptide_set}/batches/batch_{batch_id}.csv",
        tool=tool,
        peptide_set=wildcards.peptide_set,
        batch_id=batch_ids,
    )
    return [ancient(report) for report in reports]

rule algpred2_batch:
    input:
        fasta=non_hemo_batch_fasta,
        mapping=non_hemo_batch_mapping,
    output:
        report=(
            "results/inmuno_check/algpred2/{peptide_set}/batches/"
            "batch_{batch_id}.csv"
        ),
    params:
        workdir="results/inmuno_check/algpred2/{peptide_set}/work/batch_{batch_id}",
        raw_report="results/inmuno_check/algpred2/{peptide_set}/work/batch_{batch_id}/raw_algpred2.csv",
    threads: 1
    resources:
        mem_mb=get_mem_mb
    conda:
        "../envs/inmuno_check/algpred2.yml"
    shell:
        r"""
        rm -rf {params.workdir}
        mkdir -p {params.workdir}
        output_report="$(realpath -m {output.report})"
        workdir="$(realpath -m {params.workdir})"
        mkdir -p "$(dirname "$output_report")"

        input_fasta="$(realpath {input.fasta})"

        cd "$workdir"
        algpred2 -i "$input_fasta" -o raw_algpred2.csv -d 2
        cd -

        python code/hem_check/annotate_indexed_report.py \
            --raw-report {params.raw_report} \
            --mapping {input.mapping} \
            --output "$output_report" \
            --prefix algpred2
        """

rule merge_algpred2_batches:
    input:
        lambda wildcards: inmuno_check_batch_reports(wildcards, "algpred2")
    output:
        report=(
            "results/inmuno_check/algpred2/{peptide_set}/"
            "clusters_{peptide_set}_rep_seq_algpred2.csv"
        ),
    conda:
        "../envs/tox_check/toxinpred3_captp.yml"
    shell:
        r"""
        python code/merge_csv_reports.py \
            --inputs {input} \
            --output {output.report}
        """

rule download_allergenai_model:
    output:
        model="code/inmuno_check/finalmodel.h5"
    shell:
        "wget -qO {output.model} https://compbio.uth.edu/AllergenAI/finalmodel.h5"

rule allergenai_batch:
    input:
        fasta=non_hemo_batch_fasta,
        mapping=non_hemo_batch_mapping,
        model="code/inmuno_check/finalmodel.h5",
    output:
        report=(
            "results/inmuno_check/allergenai/{peptide_set}/batches/"
            "batch_{batch_id}.csv"
        ),
    params:
        workdir="results/inmuno_check/allergenai/{peptide_set}/work/batch_{batch_id}",
        raw_report="results/inmuno_check/allergenai/{peptide_set}/work/batch_{batch_id}/raw_allergenai.csv",
    threads: 1
    resources:
        mem_mb=get_mem_mb
    conda:
        "../envs/inmuno_check/allergenai.yml"
    shell:
        r"""
        rm -rf {params.workdir}
        mkdir -p {params.workdir}
        
        # Copiar input fasta al workdir
        cp {input.fasta} {params.workdir}/input.fasta
        
        # Ejecutar preprocesamiento (creará carpetas y archivos en {params.workdir})
        python code/inmuno_check/AllergenAI_preprocess.py {params.workdir}/input.fasta
        
        # Ejecutar el modelo
        export ALLERGENAI_MODEL_PATH=code/inmuno_check/finalmodel.h5
        python code/inmuno_check/Run_AllergenAI.py {params.workdir}/input.txt {params.workdir}/input_idx.txt {params.raw_report}

        # Anotar resultados
        python code/hem_check/annotate_indexed_report.py \
            --raw-report {params.raw_report} \
            --mapping {input.mapping} \
            --output {output.report} \
            --prefix allergenai
        """

rule merge_allergenai_batches:
    input:
        lambda wildcards: inmuno_check_batch_reports(wildcards, "allergenai")
    output:
        report=(
            "results/inmuno_check/allergenai/{peptide_set}/"
            "clusters_{peptide_set}_rep_seq_allergenai.csv"
        ),
    conda:
        "../envs/tox_check/toxinpred3_captp.yml"
    shell:
        r"""
        python code/merge_csv_reports.py \
            --inputs {input} \
            --output {output.report}
        """

rule allertrans_batch:
    input:
        fasta=non_hemo_batch_fasta,
        mapping=non_hemo_batch_mapping,
    output:
        report=(
            "results/inmuno_check/allertrans/{peptide_set}/batches/"
            "batch_{batch_id}.csv"
        ),
    params:
        workdir="results/inmuno_check/allertrans/{peptide_set}/work/batch_{batch_id}",
        raw_report="results/inmuno_check/allertrans/{peptide_set}/work/batch_{batch_id}/raw_allertrans.csv",
    threads: 1
    resources:
        mem_mb=get_mem_mb
    conda:
        "../envs/inmuno_check/allertrans.yml"
    shell:
        r"""
        rm -rf {params.workdir}
        mkdir -p {params.workdir}
        
        # Copiar input fasta al workdir
        input_fasta="$(realpath {input.fasta})"
        raw_report="$(realpath {params.raw_report})"
        
        # Ejecutar AllerTrans
        cd code/inmuno_check/AllerTrans/src
        python run_all.py --fasta "$input_fasta" --output "$raw_report"
        cd -
        
        # Anotar resultados
        python code/hem_check/annotate_indexed_report.py \
            --raw-report {params.raw_report} \
            --mapping {input.mapping} \
            --output {output.report} \
            --prefix allertrans
        """

rule merge_allertrans_batches:
    input:
        lambda wildcards: inmuno_check_batch_reports(wildcards, "allertrans")
    output:
        report=(
            "results/inmuno_check/allertrans/{peptide_set}/"
            "clusters_{peptide_set}_rep_seq_allertrans.csv"
        ),
    conda:
        "../envs/tox_check/toxinpred3_captp.yml"
    shell:
        r"""
        python code/merge_csv_reports.py \
            --inputs {input} \
            --output {output.report}
        """


def get_allergenai_input(wildcards):
    if wildcards.peptide_set == "10_25":
        return []
    return ancient(
        f"results/inmuno_check/allergenai/{wildcards.peptide_set}/"
        f"clusters_{wildcards.peptide_set}_rep_seq_allergenai.csv"
    )

rule build_inmuno_summary:
    input:
        script="code/inmuno_check/build_inmuno_summary.py",
        fasta=(
            f"{NON_HEMO_FASTA_ROOT}/{{peptide_set}}/"
            "clusters_{peptide_set}_rep_seq_non_hemo.fasta"
        ),
        algpred2=ancient(
            "results/inmuno_check/algpred2/{peptide_set}/"
            "clusters_{peptide_set}_rep_seq_algpred2.csv"
        ),
        allergenai=get_allergenai_input,
        allertrans=ancient(
            "results/inmuno_check/allertrans/{peptide_set}/"
            "clusters_{peptide_set}_rep_seq_allertrans.csv"
        ),
    output:
        summary=(
            "results/inmuno_check/inmuno_summary/{peptide_set}/"
            "clusters_{peptide_set}_inmuno_summary.csv"
        ),
    conda:
        "../envs/tox_check/toxinpred3_captp.yml"
    shell:
        r"""
        python {input.script} \
            --fasta {input.fasta} \
            --algpred2 {input.algpred2} \
            --allertrans {input.allertrans} \
            $(if [ -n "{input.allergenai}" ]; then echo "--allergenai {input.allergenai}"; fi) \
            --output-csv {output.summary}
        """

rule filter_post_checks_fasta:
    input:
        script="code/inmuno_check/filter_fasta_by_inmuno_summary.py",
        fasta=(
            f"{NON_HEMO_FASTA_ROOT}/{{peptide_set}}/"
            "clusters_{peptide_set}_rep_seq_non_hemo.fasta"
        ),
        summary=ancient(
            "results/inmuno_check/inmuno_summary/{peptide_set}/"
            "clusters_{peptide_set}_inmuno_summary.csv"
        ),
    output:
        fasta=(
            f"{POST_CHECKS_FASTA_ROOT}/{{peptide_set}}/"
            "clusters_{peptide_set}_rep_seq_post_checks.fasta"
        ),
        stats=(
            f"{POST_CHECKS_FASTA_ROOT}/{{peptide_set}}/"
            "clusters_{peptide_set}_rep_seq_post_checks.stats"
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
