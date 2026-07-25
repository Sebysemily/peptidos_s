import argparse
import csv
from collections import defaultdict
from pathlib import Path

def count_fasta(path):
    count = 0
    with open(path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                count += 1
    return count

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--initial", nargs="+", required=True)
    parser.add_argument("--non-toxic", nargs="+", required=True)
    parser.add_argument("--non-hemo", nargs="+", required=True)
    parser.add_argument("--post-checks", nargs="+", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    # Map file paths to their parent's parent dir name to get peptide_set, or just assume order matches.
    # Actually, the file names contain the peptide set: clusters_{peptide_set}_rep_seq...
    def extract_set(path):
        return Path(path).parent.name
    
    data = []
    
    for fasta in args.initial:
        pset = extract_set(fasta)
        data.append({"peptide_set": pset, "step": "1_initial", "count": count_fasta(fasta)})
        
    for fasta in args.non_toxic:
        pset = extract_set(fasta)
        data.append({"peptide_set": pset, "step": "2_non_toxic", "count": count_fasta(fasta)})
        
    for fasta in args.non_hemo:
        pset = extract_set(fasta)
        data.append({"peptide_set": pset, "step": "3_non_hemo", "count": count_fasta(fasta)})
        
    for fasta in args.post_checks:
        pset = extract_set(fasta)
        data.append({"peptide_set": pset, "step": "4_post_checks", "count": count_fasta(fasta)})
        
    data.sort(key=lambda x: (x["peptide_set"], x["step"]))
    
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["peptide_set", "step", "count"])
        writer.writeheader()
        writer.writerows(data)

if __name__ == "__main__":
    main()
