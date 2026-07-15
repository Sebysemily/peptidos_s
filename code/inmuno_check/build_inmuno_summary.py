#!/usr/bin/env python3
import argparse
import csv
from collections import defaultdict, deque
from pathlib import Path


def read_fasta_records(path):
    records = []
    header = None
    seq_lines = []
    with open(path, "r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    sequence = "".join(part.strip() for part in seq_lines)
                    if sequence:
                        records.append((header[1:], sequence))
                header = line
                seq_lines = []
            else:
                seq_lines.append(line)
    if header is not None:
        sequence = "".join(part.strip() for part in seq_lines)
        if sequence:
            records.append((header[1:], sequence))
    return records


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


def is_allergen(prediction):
    if not prediction:
        return False
    pred_lower = prediction.lower()
    if "non-allergen" in pred_lower or "non allergen" in pred_lower:
        return False
    if "allergen" in pred_lower:
        return True
    return False


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--algpred2", required=True)
    parser.add_argument("--allergenai", required=False, default=None)
    parser.add_argument("--allertrans", required=True)
    parser.add_argument("--output-csv", required=True)
    args = parser.parse_args()

    fasta_records = read_fasta_records(args.fasta)
    algpred2_by_sequence = rows_by_sequence(read_csv(args.algpred2), "sequence")
    
    allergenai_by_sequence = {}
    if args.allergenai:
        allergenai_by_sequence = rows_by_sequence(read_csv(args.allergenai), "sequence")
        
    allertrans_by_sequence = rows_by_sequence(read_csv(args.allertrans), "sequence")

    output_rows = []

    for index, (header, sequence) in enumerate(fasta_records, start=1):
        alg_row = pop_sequence_row(algpred2_by_sequence, sequence)
        ai_row = pop_sequence_row(allergenai_by_sequence, sequence)
        trans_row = pop_sequence_row(allertrans_by_sequence, sequence)

        alg_pred = alg_row.get("algpred2_prediction", "") if alg_row else ""
        ai_pred = ai_row.get("allergenai_prediction", "") if ai_row else ""
        trans_pred = trans_row.get("allertrans_prediction", "") if trans_row else ""

        votes = 0
        if is_allergen(alg_pred): votes += 1
        if is_allergen(ai_pred): votes += 1
        if is_allergen(trans_pred): votes += 1

        passes_filter = votes < 2

        row = {
            "peptide_index": index,
            "peptide_id": header,
            "sequence": sequence,
            "length": len(sequence),
            "algpred2_prediction": alg_pred,
            "allergenai_prediction": ai_pred,
            "allertrans_prediction": trans_pred,
            "allergen_votes": votes,
            "inmuno_filter_pass": str(passes_filter).lower(),
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
