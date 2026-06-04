#!/usr/bin/env python3
import argparse
import csv
from collections import defaultdict, deque
from pathlib import Path


def read_fasta_records(path):
    records = []
    header = None
    seq_lines = []

    with open(path, "r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    sequence = "".join(part.strip() for part in seq_lines)
                    if sequence:
                        records.append((header[1:], sequence))
                header = line
                seq_lines = []
            else:
                seq_lines.append(line)

    if header is not None:
        sequence = "".join(part.strip() for part in seq_lines)
        if sequence:
            records.append((header[1:], sequence))

    return records


def read_csv(path):
    with open(path, "r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def rows_by_sequence(rows, sequence_column):
    grouped = defaultdict(deque)
    for row in rows:
        sequence = row.get(sequence_column, "").strip()
        if sequence:
            grouped[sequence].append(row)
    return grouped


def pop_sequence_row(grouped, sequence):
    if sequence not in grouped or not grouped[sequence]:
        return None
    return grouped[sequence].popleft()


def main():
    parser = argparse.ArgumentParser(
        description="Build a per-peptide hemolytic summary CSV."
    )
    parser.add_argument("--fasta", required=True, help="Representative FASTA")
    parser.add_argument("--macrel", required=True, help="Merged Macrel CSV")
    parser.add_argument("--hemopi2-classification", required=True, help="Merged HemoPI-2 Classification CSV")
    parser.add_argument("--hemopi2-regression", required=True, help="Merged HemoPI-2 Regression CSV")
    parser.add_argument("--hepad-hmp1", required=True, help="Merged HEPAD Hmp1 CSV")
    parser.add_argument("--hepad-hmpm", required=True, help="Merged HEPAD Hmpm CSV")
    parser.add_argument("--output-csv", required=True, help="Output summary CSV")
    args = parser.parse_args()

    fasta_records = read_fasta_records(args.fasta)

    macrel_by_seq = rows_by_sequence(read_csv(args.macrel), "sequence")
    hemopi2_cls_by_seq = rows_by_sequence(read_csv(args.hemopi2_classification), "sequence")
    hemopi2_reg_by_seq = rows_by_sequence(read_csv(args.hemopi2_regression), "sequence")
    hepad_hmp1_by_seq = rows_by_sequence(read_csv(args.hepad_hmp1), "sequence")
    hepad_hmpm_by_seq = rows_by_sequence(read_csv(args.hepad_hmpm), "sequence")

    output_rows = []

    for index, (header, sequence) in enumerate(fasta_records, start=1):
        macrel_row = pop_sequence_row(macrel_by_seq, sequence) or {}
        hemopi2_cls_row = pop_sequence_row(hemopi2_cls_by_seq, sequence) or {}
        hemopi2_reg_row = pop_sequence_row(hemopi2_reg_by_seq, sequence) or {}
        hepad_hmp1_row = pop_sequence_row(hepad_hmp1_by_seq, sequence) or {}
        hepad_hmpm_row = pop_sequence_row(hepad_hmpm_by_seq, sequence) or {}

        row = {
            "peptide_index": index,
            "peptide_id": header,
            "sequence": sequence,
            "length": len(sequence),
            "macrel_hemolytic_prediction": macrel_row.get("macrel_hemolytic_prediction", ""),
            "hemopi2_classification_prediction": hemopi2_cls_row.get("hemopi2_classification_prediction", ""),
            "hemopi2_regression_prediction": hemopi2_reg_row.get("hemopi2_regression_prediction", ""),
            "hepad_Hmp1_rbfsvm_binary": hepad_hmp1_row.get("hepad_Hmp1_rbfsvm_binary", ""),
            "hepad_Hmpm_lightgbm_binary": hepad_hmpm_row.get("hepad_Hmpm_lightgbm_binary", ""),
            "hepad_Hmpm_gbc_binary": hepad_hmpm_row.get("hepad_Hmpm_gbc_binary", ""),
        }
        output_rows.append(row)

    output_csv = Path(args.output_csv)
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with output_csv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output_rows[0].keys()))
        writer.writeheader()
        writer.writerows(output_rows)


if __name__ == "__main__":
    main()
