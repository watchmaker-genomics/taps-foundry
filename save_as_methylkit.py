#!/usr/bin/env python3

import pandas as pd
import argparse
import os
import sys


def generate_output_filename(input_path: str) -> str:
    base, ext = os.path.splitext(input_path)
    return f"{base}_methylkit.tsv"


def save_as_methylkit(df: pd.DataFrame, output_filename: str) -> None:
    """Convert a methylation DataFrame to methylKit format and save it."""
    try:
        df["chrBase"] = df["#chr"].astype(str) + ":" + df["start"].astype(str)
        df["chr"] = df["#chr"]
        df["base"] = df["start"]
        df["strand"] = df.get("strand", ".")
        df["strand"] = df["strand"].apply(
            lambda x: "R" if x == "-" else ("F" if x == "+" else ".")
        )
        df["coverage"] = df["coverage"]
        df["freqC"] = 100 * df["unmod"] / (df["mod"] + df["unmod"])
        df["freqT"] = 100 * df["mod"] / (df["mod"] + df["unmod"])
        df = df[["chrBase", "chr", "base", "strand", "coverage", "freqC", "freqT"]]
        df.to_csv(output_filename, sep="\t", index=False)
    except Exception as e:
        print(f"Error while converting to methylKit format: {e}", file=sys.stderr)
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description="Convert Rastair methylation output to methylKit format."
    )
    parser.add_argument("input_file", help="Path to the Rastair TSV file.")
    args = parser.parse_args()

    try:
        df = pd.read_csv(args.input_file, sep="\t")
    except Exception as e:
        print(f"Failed to read input file: {e}", file=sys.stderr)
        sys.exit(1)

    output_filename = generate_output_filename(args.input_file)
    save_as_methylkit(df, output_filename)
    print(f"Output saved to: {output_filename}")


if __name__ == "__main__":
    main()
