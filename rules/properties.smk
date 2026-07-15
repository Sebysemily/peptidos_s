PROCESS_GROUPS_TOGETHER = config.get("properties", {}).get("process_groups_together", True)

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
