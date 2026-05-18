#!/usr/bin/env python3
import argparse
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(
        description="Merge CSV reports, keeping only the first header."
    )
    parser.add_argument("--inputs", nargs="+", required=True, help="Input CSV files")
    parser.add_argument("--output", required=True, help="Merged output CSV")
    args = parser.parse_args()

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)

    with output.open("w", encoding="utf-8") as out_handle:
        for index, input_path in enumerate(args.inputs):
            path = Path(input_path)
            if not path.exists():
                raise SystemExit(f"Missing input report: {path}")

            with path.open("r", encoding="utf-8") as in_handle:
                for line_number, line in enumerate(in_handle):
                    if index > 0 and line_number == 0:
                        continue
                    out_handle.write(line)


if __name__ == "__main__":
    main()
