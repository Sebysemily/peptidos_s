#!/usr/bin/env python3
import argparse
import subprocess
import pandas as pd
from pathlib import Path
import re
import os

def parse_anticp2_output(raw_output_path, model_name):
    records = []
    if not Path(raw_output_path).exists():
        return pd.DataFrame()
        
    with open(raw_output_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            # Remove leading '>' from Sequence_ID if present
            if line.startswith(">"):
                line = line[1:]

            # Split by comma as the output is CSV format
            parts = line.split(',')
            if len(parts) >= 4:
                seq_id = parts[0].strip()
                seq = parts[1].strip()
                
                # "score es la ultima y prediction la penultima"
                prediction = parts[-2].strip()
                score = parts[-1].strip()
                
                records.append({
                    "indexed_id": seq_id,
                    "sequence": seq,
                    f"score_model{model_name}": score,
                    f"prediction_model{model_name}": prediction,
                })

    return pd.DataFrame(records)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--raw-output", required=True)
    parser.add_argument("--csv-output", required=True)
    parser.add_argument("--mapping", required=False)
    args = parser.parse_args()

    # Paths for temporary raw outputs
    raw_path = Path(args.raw_output)
    raw_m1 = raw_path.parent / f"m1_{raw_path.name}"
    raw_m2 = raw_path.parent / f"m2_{raw_path.name}"

    # Run Model 1 with -d 2 to get all peptides (not just anticancer)
    subprocess.run(["anticp2", "-i", args.fasta, "-o", str(raw_m1), "-m", "1", "-d", "2"], check=True)
    # Run Model 2 with -d 2
    subprocess.run(["anticp2", "-i", args.fasta, "-o", str(raw_m2), "-m", "2", "-d", "2"], check=True)

    df1 = parse_anticp2_output(raw_m1, "1")
    df2 = parse_anticp2_output(raw_m2, "2")
    
    if df1.empty and df2.empty:
        # Create empty df with expected columns
        df = pd.DataFrame(columns=["indexed_id", "sequence", "prediction_model1", "score_model1", "prediction_model2", "score_model2"])
    elif df1.empty:
        df = df2
    elif df2.empty:
        df = df1
    else:
        # Merge both outputs on indexed_id and sequence
        df = pd.merge(df1, df2, on=["indexed_id", "sequence"], how="outer")

    # Restore original IDs using mapping
    if args.mapping and Path(args.mapping).exists():
        mapping_df = pd.read_csv(args.mapping)
        if not df.empty:
            df = df.merge(
                mapping_df[["indexed_id", "peptide_id"]], on="indexed_id", how="left"
            )
            df = df.drop(columns=["indexed_id"])
            # Reorder columns
            cols = ["peptide_id"] + [c for c in df.columns if c != "peptide_id"]
            df = df[cols]
        else:
            df = pd.DataFrame(columns=["peptide_id", "sequence", "prediction_model1", "score_model1", "prediction_model2", "score_model2"])
    else:
        if not df.empty:
            df = df.rename(columns={"indexed_id": "peptide_id"})
            cols = ["peptide_id"] + [c for c in df.columns if c != "peptide_id"]
            df = df[cols]

    df.to_csv(args.csv_output, index=False)
    
    # Cleanup raw files
    if raw_m1.exists():
        raw_m1.unlink()
    if raw_m2.exists():
        raw_m2.unlink()


if __name__ == "__main__":
    main()
