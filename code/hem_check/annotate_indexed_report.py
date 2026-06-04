#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path


ID_COLUMN_HINTS = ("id", "name", "header", "accession", "seq")


def read_mapping(path):
    rows = []
    by_id = {}
    by_sequence = {}
    with open(path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        id_field = "indexed_id" if "indexed_id" in (reader.fieldnames or []) else "tox_index"
        for row in reader:
            indexed_id = row[id_field]
            row["indexed_id"] = indexed_id
            rows.append(row)
            by_id[indexed_id] = row
            by_sequence.setdefault(row["sequence"], []).append(row)
    return rows, by_id, by_sequence


def find_mapping(row, by_id, by_sequence):
    hinted_columns = [
        column for column in row
        if any(hint in column.lower().replace("_", " ") for hint in ID_COLUMN_HINTS)
    ]
    for column in hinted_columns + list(row):
        value = (row.get(column) or "").strip()
        if value in by_id:
            return by_id[value]
        if value.startswith(">") and value[1:] in by_id:
            return by_id[value[1:]]

    for value in row.values():
        matches = by_sequence.get((value or "").strip())
        if matches and len(matches) == 1:
            return matches[0]

    return None


def prefixed_fieldnames(fieldnames, prefix):
    clean = []
    for field in fieldnames:
        if field in {"indexed_id", "peptide_id", "sequence", "length"}:
            clean.append(f"{prefix}_{field}")
        else:
            clean.append(f"{prefix}_{field.strip().lower().replace(' ', '_')}")
    return clean


def main():
    parser = argparse.ArgumentParser(
        description="Annotate a report produced from indexed FASTA input."
    )
    parser.add_argument("--raw-report", required=True, help="Raw tool report CSV")
    parser.add_argument("--mapping", required=True, help="Indexed FASTA mapping CSV")
    parser.add_argument("--output", required=True, help="Annotated output CSV")
    parser.add_argument("--prefix", required=True, help="Prefix for raw report columns")
    args = parser.parse_args()

    _, by_id, by_sequence = read_mapping(args.mapping)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)

    with open(args.raw_report, "r", encoding="utf-8", newline="") as raw_handle:
        reader = csv.DictReader(raw_handle)
        raw_fields = reader.fieldnames or []
        renamed_fields = prefixed_fieldnames(raw_fields, args.prefix)
        fieldnames = ["indexed_id", "peptide_id", "sequence", "length"] + renamed_fields

        with output.open("w", encoding="utf-8", newline="") as output_handle:
            writer = csv.DictWriter(output_handle, fieldnames=fieldnames)
            writer.writeheader()

            for row_number, row in enumerate(reader, start=2):
                mapped = find_mapping(row, by_id, by_sequence)
                if mapped is None:
                    raise SystemExit(
                        f"Could not map row {row_number} from {args.raw_report} "
                        "to an indexed FASTA record"
                    )

                output_row = {
                    "indexed_id": mapped["indexed_id"],
                    "peptide_id": mapped["peptide_id"],
                    "sequence": mapped["sequence"],
                    "length": mapped["length"],
                }
                for raw_field, renamed_field in zip(raw_fields, renamed_fields):
                    output_row[renamed_field] = row.get(raw_field, "")
                writer.writerow(output_row)


if __name__ == "__main__":
    main()
