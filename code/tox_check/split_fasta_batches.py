#!/usr/bin/env python3
import argparse
import os
import math
import shutil
from pathlib import Path


def read_fasta_records(path):
    records = []
    header = None
    seq_lines = []
    skipped_empty = 0

    with open(path, "r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    if sequence_length(seq_lines) > 0:
                        records.append((header, seq_lines))
                    else:
                        skipped_empty += 1
                header = line
                seq_lines = []
            else:
                seq_lines.append(line)

    if header is not None:
        if sequence_length(seq_lines) > 0:
            records.append((header, seq_lines))
        else:
            skipped_empty += 1

    return records, skipped_empty


def sequence_length(seq_lines):
    return sum(len(line.strip()) for line in seq_lines)


def write_records(path, records):
    with open(path, "w", encoding="utf-8") as handle:
        for header, seq_lines in records:
            handle.write(f"{header}\n")
            for line in seq_lines:
                handle.write(f"{line}\n")


def batch_bounds(n_records, n_batches):
    base = n_records // n_batches
    remainder = n_records % n_batches
    start = 0
    for batch_index in range(n_batches):
        size = base + (1 if batch_index < remainder else 0)
        end = start + size
        yield start, end
        start = end


def fixed_size_bounds(n_records, max_records_per_batch):
    for start in range(0, n_records, max_records_per_batch):
        yield start, min(start + max_records_per_batch, n_records)


def main():
    parser = argparse.ArgumentParser(
        description="Split a FASTA file into contiguous batches."
    )
    parser.add_argument("--input", required=True, help="Input FASTA file")
    parser.add_argument("--outdir", required=True, help="Output batch directory")
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument(
        "--n-batches", type=int, help="Number of near-equal FASTA batches"
    )
    mode.add_argument(
        "--max-records-per-batch",
        type=int,
        help="Maximum number of FASTA records per batch",
    )
    args = parser.parse_args()

    records, skipped_empty = read_fasta_records(args.input)
    if not records:
        raise SystemExit(f"No non-empty FASTA records found in {args.input}")

    if args.n_batches is not None:
        if args.n_batches < 1:
            raise SystemExit("--n-batches must be at least 1")
        if args.n_batches > len(records):
            raise SystemExit(
                f"Requested {args.n_batches} batches for {len(records)} records"
            )
        bounds = list(batch_bounds(len(records), args.n_batches))
        n_batches = args.n_batches
    else:
        if args.max_records_per_batch < 1:
            raise SystemExit("--max-records-per-batch must be at least 1")
        bounds = list(fixed_size_bounds(len(records), args.max_records_per_batch))
        n_batches = math.ceil(len(records) / args.max_records_per_batch)

    outdir = Path(args.outdir)
    if outdir.exists():
        shutil.rmtree(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    for batch_number, (start, end) in enumerate(bounds, start=1):
        write_records(outdir / f"batch_{batch_number}.fasta", records[start:end])

    complete = outdir / ".complete"
    complete.write_text(
        f"input={os.path.abspath(args.input)}\n"
        f"records={len(records)}\n"
        f"skipped_empty={skipped_empty}\n"
        f"batches={n_batches}\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
