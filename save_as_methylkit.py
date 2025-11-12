#!/usr/bin/env python3

import pandas as pd
import argparse
import os
import sys


def generate_output_filename(input_path: str) -> str:
    """
    Given an input file path, return a new filename with '_methylkit.tsv' appended.
    """
    base, _ = os.path.splitext(input_path)
    return f"{base}_methylkit.tsv"


def save_as_methylkit(df: pd.DataFrame, output_filename: str) -> None:
    """
    Convert a Rastair mods file to methylKit format and save it.

    Required input columns:
      #chr, start, mod, unmod, coverage, optional 'strand'
    """
    try:
        df["chrBase"] = df["#chr"].astype(str) + ":" + df["start"].astype(str)
        df["chr"] = df["#chr"]
        df["base"] = df["start"]

        # Convert strand to methylKit notation (R/F/.)
        df["strand"] = df.get("strand", ".")
        df["strand"] = df["strand"].apply(
            lambda x: "R" if x == "-" else ("F" if x == "+" else ".")
        )

        total = df["mod"] + df["unmod"]
        df["freqC"] = (100 * df["mod"] / total).where(total != 0, 0)
        df["freqT"] = (100 * df["unmod"] / total).where(total != 0, 0)

        df = df[["chrBase", "chr", "base", "strand", "coverage", "freqC", "freqT"]]
        df.to_csv(output_filename, sep="\t", index=False)

    except Exception as e:
        print(f"Error while converting to methylKit format: {e}", file=sys.stderr)
        sys.exit(1)


def main() -> None:
    """
    Parse arguments, load input TSV, convert, and save methylKit-formatted output.
    """
    parser = argparse.ArgumentParser(
        description="Convert Rastair methylation output to methylKit format."
    )

    parser.add_argument(
        "-i", "--input_file",
        type=str,
        required=True,
        help="Path to the Rastair TSV file to convert."
    )

    parser.add_argument(
        "-o", "--output_file",
        type=str,
        default=None,
        help="Optional output file path. If omitted, a filename will be generated automatically."
    )

    args = parser.parse_args()

    try:
        df = pd.read_csv(args.input_file, sep="\t")
    except Exception as e:
        print(f"Failed to read input file: {e}", file=sys.stderr)
        sys.exit(1)

    output_filename = args.output_file or generate_output_filename(args.input_file)

    save_as_methylkit(df, output_filename)

    print(f"Output saved to: {output_filename}")

if __name__ == "__main__":
    main()
