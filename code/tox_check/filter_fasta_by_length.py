#!/usr/bin/env python3
import argparse
from pathlib import Path


def read_fasta(path):
    header = None
    seq_lines = []

    with open(path, "r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, seq_lines
                header = line
                seq_lines = []
            else:
                seq_lines.append(line)

    if header is not None:
        yield header, seq_lines


def sequence_length(seq_lines):
    return sum(len(line.strip()) for line in seq_lines)


def write_record(handle, header, seq_lines):
    handle.write(f"{header}\n")
    for line in seq_lines:
        handle.write(f"{line}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Filter FASTA records by sequence length."
    )
    parser.add_argument("--input", required=True, help="Input FASTA")
    parser.add_argument("--output", required=True, help="Output FASTA")
    parser.add_argument("--stats", required=True, help="Output stats file")
    parser.add_argument("--min-length", type=int, default=1)
    parser.add_argument("--max-length", type=int, required=True)
    args = parser.parse_args()

    output = Path(args.output)
    stats = Path(args.stats)
    output.parent.mkdir(parents=True, exist_ok=True)
    stats.parent.mkdir(parents=True, exist_ok=True)

    total = 0
    kept = 0
    removed = 0
    removed_too_short = 0
    removed_too_long = 0

    with output.open("w", encoding="utf-8") as out_handle:
        for header, seq_lines in read_fasta(args.input):
            total += 1
            length = sequence_length(seq_lines)
            if args.min_length <= length <= args.max_length:
                kept += 1
                write_record(out_handle, header, seq_lines)
            else:
                removed += 1
                if length < args.min_length:
                    removed_too_short += 1
                else:
                    removed_too_long += 1

    stats.write_text(
        f"input={Path(args.input).resolve()}\n"
        f"min_length={args.min_length}\n"
        f"max_length={args.max_length}\n"
        f"total={total}\n"
        f"kept={kept}\n"
        f"removed={removed}\n"
        f"removed_too_short={removed_too_short}\n"
        f"removed_too_long={removed_too_long}\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
