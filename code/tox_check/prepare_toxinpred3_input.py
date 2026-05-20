#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path


def read_fasta(path):
    header = None
    seq_lines = []

    with open(path, "r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    sequence = "".join(part.strip() for part in seq_lines)
                    if sequence:
                        yield header[1:], sequence
                header = line
                seq_lines = []
            else:
                seq_lines.append(line)

    if header is not None:
        sequence = "".join(part.strip() for part in seq_lines)
        if sequence:
            yield header[1:], sequence


def wrap_sequence(sequence, width=80):
    for start in range(0, len(sequence), width):
        yield sequence[start : start + width]


def main():
    parser = argparse.ArgumentParser(
        description="Create indexed FASTA and mapping for ToxinPred3."
    )
    parser.add_argument("--input", required=True, help="Input FASTA")
    parser.add_argument("--output-fasta", required=True, help="Indexed FASTA")
    parser.add_argument("--mapping", required=True, help="ID mapping CSV")
    parser.add_argument("--prefix", required=True, help="Indexed ID prefix")
    args = parser.parse_args()

    output_fasta = Path(args.output_fasta)
    mapping = Path(args.mapping)
    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    mapping.parent.mkdir(parents=True, exist_ok=True)

    with output_fasta.open("w", encoding="utf-8") as fasta_handle, mapping.open(
        "w", encoding="utf-8", newline=""
    ) as mapping_handle:
        writer = csv.DictWriter(
            mapping_handle,
            fieldnames=["tox_index", "peptide_id", "sequence", "length"],
        )
        writer.writeheader()

        for index, (peptide_id, sequence) in enumerate(read_fasta(args.input), start=1):
            tox_index = f"{args.prefix}_{index}"
            fasta_handle.write(f">{tox_index}\n")
            for line in wrap_sequence(sequence):
                fasta_handle.write(f"{line}\n")
            writer.writerow(
                {
                    "tox_index": tox_index,
                    "peptide_id": peptide_id,
                    "sequence": sequence,
                    "length": len(sequence),
                }
            )


if __name__ == "__main__":
    main()
