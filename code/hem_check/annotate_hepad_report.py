#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path


MODELS = ("rbfsvm", "lightgbm", "gbc")
BASE_FIELDS = ("indexed_id", "peptide_id", "sequence", "length")


def read_mapping(path):
    mapping = {}
    with open(path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        id_field = "indexed_id" if "indexed_id" in (reader.fieldnames or []) else "tox_index"
        for row in reader:
            mapping[row[id_field]] = row
    return mapping


def read_prediction_table(path):
    rows = {}
    with open(path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader, None)
        if header is None:
            raise SystemExit(f"Empty prediction table: {path}")

        if header[0] == "" or header[0].lower() in {"index", "id", "name"}:
            model_names = header[1:]
            id_index = 0
        else:
            model_names = header[1:]
            id_index = 0

        for line in reader:
            if not line:
                continue
            record_id = line[id_index].strip()
            values = line[1:]
            rows[record_id] = dict(zip(model_names, values))
    return rows, model_names


def output_fields(prefix):
    fields = list(BASE_FIELDS)
    for model in MODELS:
        fields.append(f"{prefix}_{model}_binary")
    for model in MODELS:
        fields.append(f"{prefix}_{model}_probability")
    return fields


def main():
    parser = argparse.ArgumentParser(
        description="Annotate HEPAD prediction tables with original peptide IDs."
    )
    parser.add_argument("--binary-report", required=True, help="binary_vector.csv")
    parser.add_argument(
        "--probability-report", required=True, help="probability_vector.csv"
    )
    parser.add_argument("--mapping", required=True, help="Indexed FASTA mapping CSV")
    parser.add_argument("--output", required=True, help="Annotated output CSV")
    parser.add_argument(
        "--prefix",
        required=True,
        help="Prefix for HEPAD output columns (for example hepad_Hmp1)",
    )
    args = parser.parse_args()

    mapping = read_mapping(args.mapping)
    binary_rows, binary_models = read_prediction_table(args.binary_report)
    probability_rows, probability_models = read_prediction_table(args.probability_report)

    if list(binary_models) != list(MODELS):
        raise SystemExit(
            f"Unexpected binary model columns in {args.binary_report}: {binary_models}"
        )
    if list(probability_models) != list(MODELS):
        raise SystemExit(
            "Unexpected probability model columns in "
            f"{args.probability_report}: {probability_models}"
        )

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fields = output_fields(args.prefix)

    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()

        for indexed_id in mapping:
            mapped = mapping[indexed_id]
            binary = binary_rows.get(indexed_id)
            probability = probability_rows.get(indexed_id)
            if binary is None or probability is None:
                raise SystemExit(
                    f"HEPAD prediction missing for indexed_id {indexed_id!r}"
                )

            row = {
                "indexed_id": mapped["indexed_id"],
                "peptide_id": mapped["peptide_id"],
                "sequence": mapped["sequence"],
                "length": mapped["length"],
            }
            for model in MODELS:
                row[f"{args.prefix}_{model}_binary"] = binary.get(model, "")
            for model in MODELS:
                row[f"{args.prefix}_{model}_probability"] = probability.get(model, "")
            writer.writerow(row)


if __name__ == "__main__":
    main()
