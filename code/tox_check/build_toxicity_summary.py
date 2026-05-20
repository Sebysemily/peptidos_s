#!/usr/bin/env python3
import argparse
import csv
from collections import defaultdict, deque
from pathlib import Path


TOXTELLER_COLUMNS = [
    "LR prediction",
    "SVM prediction",
    "RF prediction",
    "XGBoost prediction",
]


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
                    sequence = "".join(part.strip() for part in seq_lines)
                    if sequence:
                        records.append((header[1:], sequence))
                    else:
                        skipped_empty += 1
                header = line
                seq_lines = []
            else:
                seq_lines.append(line)

    if header is not None:
        sequence = "".join(part.strip() for part in seq_lines)
        if sequence:
            records.append((header[1:], sequence))
        else:
            skipped_empty += 1

    return records, skipped_empty


def read_csv(path):
    with open(path, "r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def rows_by_sequence(rows, sequence_column):
    grouped = defaultdict(deque)
    for row in rows:
        sequence = row.get(sequence_column, "").strip()
        if sequence:
            grouped[sequence].append(row)
    return grouped


def pop_sequence_row(grouped, sequence):
    if sequence not in grouped or not grouped[sequence]:
        return None
    return grouped[sequence].popleft()


def normalized_prediction(value):
    return (value or "").strip().lower().replace("_", "-")


def is_toxinpred3_below_threshold(row, threshold=0.38):
    if row is None:
        return False
    try:
        return float(row.get("toxinpred3_hybrid_score", "")) < threshold
    except ValueError:
        return False


def is_toxteller_non_toxic(row):
    if row is None:
        return False
    predictions = []
    for column in TOXTELLER_COLUMNS:
        value = (row.get(column) or "").strip()
        if value == "":
            return False
        predictions.append(value)
    return all(value in {"0", "0.0"} for value in predictions)


def is_captp_non_toxic(row):
    if row is None:
        return None
    prediction = normalized_prediction(row.get("Prediction"))
    return "toxic" not in prediction


def main():
    parser = argparse.ArgumentParser(
        description="Build a per-peptide toxicity summary CSV."
    )
    parser.add_argument("--fasta", required=True, help="Representative FASTA")
    parser.add_argument("--toxinpred3", required=True, help="Merged ToxinPred3 CSV")
    parser.add_argument("--toxteller", required=True, help="Merged ToxTeller CSV")
    parser.add_argument("--captp", required=True, help="Merged CAPTP CSV")
    parser.add_argument("--output-csv", required=True, help="Output summary CSV")
    args = parser.parse_args()

    fasta_records, _ = read_fasta_records(args.fasta)
    toxinpred3_rows = read_csv(args.toxinpred3)
    if toxinpred3_rows and "sequence" not in toxinpred3_rows[0]:
        raise SystemExit(
            f"{args.toxinpred3} is an old unannotated ToxinPred3 report. "
            "Re-run the toxinpred3_batch and merge_toxinpred3_batches rules so "
            "the report includes peptide_id and sequence columns."
        )

    toxinpred3_by_sequence = rows_by_sequence(toxinpred3_rows, "sequence")
    toxteller_by_sequence = rows_by_sequence(read_csv(args.toxteller), "Sequence")
    captp_by_sequence = rows_by_sequence(read_csv(args.captp), "Sequences")

    output_rows = []

    for index, (header, sequence) in enumerate(fasta_records, start=1):
        toxinpred3_row = pop_sequence_row(toxinpred3_by_sequence, sequence)
        toxteller_row = pop_sequence_row(toxteller_by_sequence, sequence)
        captp_row = pop_sequence_row(captp_by_sequence, sequence)

        toxinpred3_pass = is_toxinpred3_below_threshold(toxinpred3_row)
        toxteller_non_toxic = is_toxteller_non_toxic(toxteller_row)
        captp_non_toxic = is_captp_non_toxic(captp_row)
        captp_available = captp_non_toxic is not None

        passes_filter = (
            toxinpred3_pass
            and toxteller_non_toxic
            and (captp_non_toxic is not False)
        )

        row = {
            "peptide_index": index,
            "peptide_id": header,
            "sequence": sequence,
            "length": len(sequence),
            "toxinpred3_prediction": (
                toxinpred3_row.get("toxinpred3_prediction", "")
                if toxinpred3_row
                else ""
            ),
            "toxinpred3_ml_score": (
                toxinpred3_row.get("toxinpred3_ml_score", "")
                if toxinpred3_row
                else ""
            ),
            "toxinpred3_hybrid_score": (
                toxinpred3_row.get("toxinpred3_hybrid_score", "")
                if toxinpred3_row
                else ""
            ),
            "toxteller_lr_prediction": (
                toxteller_row.get("LR prediction", "") if toxteller_row else ""
            ),
            "toxteller_svm_prediction": (
                toxteller_row.get("SVM prediction", "") if toxteller_row else ""
            ),
            "toxteller_rf_prediction": (
                toxteller_row.get("RF prediction", "") if toxteller_row else ""
            ),
            "toxteller_xgboost_prediction": (
                toxteller_row.get("XGBoost prediction", "") if toxteller_row else ""
            ),
            "captp_available": str(captp_available).lower(),
            "captp_prediction": captp_row.get("Prediction", "") if captp_row else "",
            "captp_confidence": captp_row.get("Confidence", "") if captp_row else "",
            "toxicity_filter_pass": str(passes_filter).lower(),
        }
        output_rows.append(row)

    output_csv = Path(args.output_csv)
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with output_csv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output_rows[0].keys()))
        writer.writeheader()
        writer.writerows(output_rows)


if __name__ == "__main__":
    main()
