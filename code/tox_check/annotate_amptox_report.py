#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-report", required=True)
    parser.add_argument("--mapping", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    # Load mapping
    by_id = {}
    with open(args.mapping, "r", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            by_id[row["indexed_id"]] = row

    # Read raw and write output
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.raw_report, "r", newline="") as f_in, open(args.output, "w", newline="") as f_out:
        reader = csv.DictReader(f_in)
        
        fieldnames = ["indexed_id", "peptide_id", "sequence", "length", "amptox_rf_prediction", "amptox_svc_prediction"]
        writer = csv.DictWriter(f_out, fieldnames=fieldnames)
        writer.writeheader()

        for row in reader:
            mapped = by_id.get(row["id"])
            if mapped:
                writer.writerow({
                    "indexed_id": mapped["indexed_id"],
                    "peptide_id": mapped["peptide_id"],
                    "sequence": mapped["sequence"],
                    "length": mapped["length"],
                    "amptox_rf_prediction": row.get("amptox_rf_prediction", ""),
                    "amptox_svc_prediction": row.get("amptox_svc_prediction", "")
                })

if __name__ == "__main__":
    main()
