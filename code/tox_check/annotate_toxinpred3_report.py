#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path


def read_mapping(path):
    with open(path, "r", encoding="utf-8", newline="") as handle:
        return {row["tox_index"]: row for row in csv.DictReader(handle)}


def main():
    parser = argparse.ArgumentParser(
        description="Annotate ToxinPred3 output with original peptide IDs."
    )
    parser.add_argument("--raw-report", required=True, help="Raw ToxinPred3 CSV")
    parser.add_argument("--mapping", required=True, help="ID mapping CSV")
    parser.add_argument("--output", required=True, help="Annotated CSV")
    args = parser.parse_args()

    mapping = read_mapping(args.mapping)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)

    with open(args.raw_report, "r", encoding="utf-8", newline="") as raw_handle:
        reader = csv.DictReader(raw_handle)
        fieldnames = [
            "peptide_id",
            "sequence",
            "length",
            "toxinpred3_subject",
            "toxinpred3_ml_score",
            "toxinpred3_merci_score_pos",
            "toxinpred3_merci_score_neg",
            "toxinpred3_hybrid_score",
            "toxinpred3_prediction",
            "toxinpred3_ppv",
        ]

        with output.open("w", encoding="utf-8", newline="") as out_handle:
            writer = csv.DictWriter(out_handle, fieldnames=fieldnames)
            writer.writeheader()

            for row in reader:
                subject = row.get("Subject", "")
                if subject not in mapping:
                    raise SystemExit(
                        f"ToxinPred3 subject {subject!r} is missing from {args.mapping}"
                    )
                mapped = mapping[subject]
                writer.writerow(
                    {
                        "peptide_id": mapped["peptide_id"],
                        "sequence": mapped["sequence"],
                        "length": mapped["length"],
                        "toxinpred3_subject": subject,
                        "toxinpred3_ml_score": row.get("ML Score", ""),
                        "toxinpred3_merci_score_pos": row.get(
                            "MERCI Score Pos", ""
                        ),
                        "toxinpred3_merci_score_neg": row.get(
                            "MERCI Score Neg", ""
                        ),
                        "toxinpred3_hybrid_score": row.get("Hybrid Score", ""),
                        "toxinpred3_prediction": row.get("Prediction", ""),
                        "toxinpred3_ppv": row.get("PPV", ""),
                    }
                )


if __name__ == "__main__":
    main()
