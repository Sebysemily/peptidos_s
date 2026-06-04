#!/usr/bin/env python3
import argparse
import csv
import gzip
from pathlib import Path


OUTPUT_FIELDS = [
    "macrel_index",
    "peptide_id",
    "sequence",
    "macrel_amp_family",
    "macrel_is_amp",
    "macrel_amp_probability",
    "macrel_hemolytic_prediction",
    "macrel_hemolytic_probability",
]


def read_mapping(path):
    mapping = {}
    with open(path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        id_field = "indexed_id" if "indexed_id" in (reader.fieldnames or []) else "tox_index"
        for row in reader:
            mapping[row[id_field]] = row
    return mapping


def prediction_file(input_dir):
    root = Path(input_dir)
    candidates = sorted(root.rglob("*.prediction"))
    candidates.extend(sorted(root.rglob("*.prediction.gz")))
    if not candidates:
        candidates = sorted(
            path
            for path in root.rglob("*")
            if path.name in {"prediction", "prediction.gz"}
        )
    if len(candidates) != 1:
        found = ", ".join(str(path) for path in candidates) or "none"
        raise SystemExit(f"Expected exactly one Macrel prediction file in {root}; found {found}")
    return candidates[0]


def split_row(line):
    line = line.rstrip("\n")
    if "\t" in line:
        return [part.strip() for part in line.split("\t")]
    return [part.strip() for part in next(csv.reader([line]))]


def is_header(parts):
    normalized = {part.strip().lower().replace("_", " ") for part in parts}
    return bool(normalized.intersection({"accession", "sequence", "seq id", "hemolytic"}))


def iter_prediction_rows(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            if line.startswith("#"):
                continue
            parts = split_row(line)
            if is_header(parts):
                continue
            if len(parts) < 7:
                raise SystemExit(
                    f"Expected at least 7 columns in {path}, got {len(parts)}: {line.rstrip()}"
                )
            yield parts[:7]


def sequence_matches(macrel_sequence, mapping_sequence):
    if macrel_sequence == mapping_sequence:
        return True
    return mapping_sequence.startswith("M") and macrel_sequence == mapping_sequence[1:]


def main():
    parser = argparse.ArgumentParser(
        description="Annotate a Macrel prediction table with original peptide IDs."
    )
    parser.add_argument("--input-dir", required=True, help="Macrel output directory")
    parser.add_argument("--mapping", required=True, help="Indexed FASTA mapping CSV")
    parser.add_argument("--output", required=True, help="Annotated output CSV")
    args = parser.parse_args()

    mapping = read_mapping(args.mapping)
    predictions = prediction_file(args.input_dir)

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)

    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_FIELDS)
        writer.writeheader()

        for (
            macrel_index,
            macrel_sequence,
            amp_family,
            is_amp,
            amp_probability,
            hemolytic,
            hemolytic_probability,
        ) in iter_prediction_rows(predictions):
            mapped = mapping.get(macrel_index)
            if mapped is None:
                raise SystemExit(f"Macrel index {macrel_index!r} is missing from {args.mapping}")

            sequence = mapped["sequence"]
            if macrel_sequence and not sequence_matches(macrel_sequence, sequence):
                raise SystemExit(
                    f"Macrel sequence mismatch for {macrel_index}: "
                    f"{macrel_sequence!r} != {sequence!r}"
                )

            writer.writerow(
                {
                    "macrel_index": macrel_index,
                    "peptide_id": mapped["peptide_id"],
                    "sequence": sequence,
                    "macrel_amp_family": amp_family,
                    "macrel_is_amp": is_amp,
                    "macrel_amp_probability": amp_probability,
                    "macrel_hemolytic_prediction": hemolytic,
                    "macrel_hemolytic_probability": hemolytic_probability,
                }
            )


if __name__ == "__main__":
    main()
