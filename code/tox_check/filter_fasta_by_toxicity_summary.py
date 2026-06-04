#!/usr/bin/env python3
import argparse
import csv
from collections import Counter
from pathlib import Path


TRUE_VALUES = {"true", "1", "yes", "y"}


def read_fasta_records(path):
    header = None
    seq_lines = []

    with open(path, "r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    sequence = "".join(part.strip() for part in seq_lines)
                    if sequence:
                        yield header[1:], seq_lines, sequence
                header = line
                seq_lines = []
            else:
                seq_lines.append(line)

    if header is not None:
        sequence = "".join(part.strip() for part in seq_lines)
        if sequence:
            yield header[1:], seq_lines, sequence


def read_passing_records(path):
    passing = Counter()
    total_rows = 0
    passing_rows = 0

    with open(path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"peptide_id", "sequence", "toxicity_filter_pass"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            missing_fields = ", ".join(sorted(missing))
            raise SystemExit(f"{path} is missing required columns: {missing_fields}")

        for row in reader:
            total_rows += 1
            passed = (row.get("toxicity_filter_pass") or "").strip().lower()
            if passed in TRUE_VALUES:
                peptide_id = row.get("peptide_id", "")
                sequence = (row.get("sequence") or "").strip()
                passing[(peptide_id, sequence)] += 1
                passing_rows += 1

    return passing, total_rows, passing_rows


def main():
    parser = argparse.ArgumentParser(
        description="Write a FASTA containing only peptides that pass toxicity filtering."
    )
    parser.add_argument("--fasta", required=True, help="Input representative FASTA")
    parser.add_argument("--summary", required=True, help="Toxicity summary CSV")
    parser.add_argument("--output-fasta", required=True, help="Filtered output FASTA")
    parser.add_argument("--stats", required=True, help="Output stats file")
    args = parser.parse_args()

    passing, total_rows, passing_rows = read_passing_records(args.summary)

    output_fasta = Path(args.output_fasta)
    stats = Path(args.stats)
    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    stats.parent.mkdir(parents=True, exist_ok=True)

    input_records = 0
    written_records = 0

    with output_fasta.open("w", encoding="utf-8") as output_handle:
        for peptide_id, seq_lines, sequence in read_fasta_records(args.fasta):
            input_records += 1
            key = (peptide_id, sequence)
            if passing[key] <= 0:
                continue

            passing[key] -= 1
            written_records += 1
            output_handle.write(f">{peptide_id}\n")
            for line in seq_lines:
                output_handle.write(f"{line}\n")

    unmatched = sum(passing.values())
    if unmatched:
        raise SystemExit(
            f"{unmatched} passing toxicity summary rows were not found in {args.fasta}"
        )

    stats.write_text(
        "\n".join(
            [
                f"input_fasta={Path(args.fasta).resolve()}",
                f"toxicity_summary={Path(args.summary).resolve()}",
                f"summary_rows={total_rows}",
                f"summary_passing_rows={passing_rows}",
                f"input_records={input_records}",
                f"written_records={written_records}",
                "",
            ]
        ),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
