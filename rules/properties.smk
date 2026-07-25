PROCESS_GROUPS_TOGETHER = config.get("properties", {}).get("process_groups_together", True)

rule count_filtered_peptides:
    input:
        initial=expand("data/curated_md-lais/mmseqs2/{peptide_set}/clusters_{peptide_set}_rep_seq.fasta", peptide_set=PEPTIDE_SETS),
        non_toxic=expand(f"{NON_TOXIC_FASTA_ROOT}/{{peptide_set}}/clusters_{{peptide_set}}_rep_seq_non_toxic.fasta", peptide_set=PEPTIDE_SETS),
        non_hemo=expand(f"{NON_HEMO_FASTA_ROOT}/{{peptide_set}}/clusters_{{peptide_set}}_rep_seq_non_hemo.fasta", peptide_set=PEPTIDE_SETS),
        post_checks=expand(f"{POST_CHECKS_FASTA_ROOT}/{{peptide_set}}/clusters_{{peptide_set}}_rep_seq_post_checks.fasta", peptide_set=PEPTIDE_SETS),
    output:
        csv="metadata/filtering.csv"
    conda:
        "../envs/properties.yml"
    shell:
        r"""
        python code/properties/count_filtered_peptides.py \
            --initial {input.initial} \
            --non-toxic {input.non_toxic} \
            --non-hemo {input.non_hemo} \
            --post-checks {input.post_checks} \
            --output {output.csv}
        """

if PROCESS_GROUPS_TOGETHER:
    rule calculate_properties_together:
        input:
            expand(
                f"{POST_CHECKS_FASTA_ROOT}/{{peptide_set}}/clusters_{{peptide_set}}_rep_seq_post_checks.fasta",
                peptide_set=PEPTIDE_SETS
            )
        output:
            "metadata/characteristics.csv"
        conda:
            "../envs/properties.yml"
        script:
            "../code/properties/calculate_properties.py"
else:
    rule calculate_properties_separate:
        input:
            f"{POST_CHECKS_FASTA_ROOT}/{{peptide_set}}/"
            "clusters_{peptide_set}_rep_seq_post_checks.fasta"
        output:
            "metadata/{peptide_set}/characteristics.csv"
        conda:
            "../envs/properties.yml"
        script:
            "../code/properties/calculate_properties.py"
