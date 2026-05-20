#!/usr/bin/env python3
import argparse
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Sanitize FASTA headers for HemoPI2's simple FASTA parser."
    )
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def main():
    args = parse_args()
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)

    with open(args.input, "r", encoding="utf-8") as in_handle, output.open(
        "w", encoding="utf-8"
    ) as out_handle:
        for raw_line in in_handle:
            line = raw_line.rstrip("\n")
            if line.startswith(">"):
                out_handle.write(">" + line[1:].replace(">", "_") + "\n")
            else:
                out_handle.write(line + "\n")


if __name__ == "__main__":
    main()
