#!/usr/bin/env python3

import io
import csv
import argparse
import os
from typing import NamedTuple
from dataclasses import dataclass


class Row(NamedTuple):
    """Represents one CpG record from the mods/summary table."""
    chr: str
    start: int
    end: int
    name: str
    beta_est: str
    strand: str
    unmod: int
    mod: int
    no_snp: int
    snp: int
    coverage: int
    genotype: str
    gt_p_score: int
    gt_conf_score: int


@dataclass
class Summary:
    """Accumulates methylation counts across all rows."""
    total_mod: int = 0
    total_unmod: int = 0
    covered_positions: int = 0

    def assimilate(self, row: Row) -> None:
        """Incorporate one row's counts into the running totals."""
        self.total_mod += row.mod
        self.total_unmod += row.unmod
        if row.coverage > 0:
            self.covered_positions += 1

    def write(self, out_fobj: io.TextIOWrapper) -> None:
        """Write a one-line summary with methylation percentage."""
        writer = csv.DictWriter(out_fobj, fieldnames=[
            "total_mod",
            "total_unmod",
            "methylation_rate",
            "covered_positions",
        ], delimiter="\t")

        denominator = self.total_mod + self.total_unmod
        methylation_rate = 0 if denominator == 0 else (self.total_mod / denominator)

        writer.writeheader()
        writer.writerow({
            "total_mod": self.total_mod,
            "total_unmod": self.total_unmod,
            "methylation_rate": round(methylation_rate, 5),
            "covered_positions": self.covered_positions,
        })


def main(input_path: str, output_path: str) -> None:
    """Generate a methylation summary from a mods file."""
    summary = Summary()

    with open(input_path, newline="") as csvfile:
        reader = csv.DictReader(csvfile, delimiter="\t")
        for row in reader:
            # Skip SNP-containing rows
            if not (row["genotype"] == "C/C" or row["genotype"] == "G/G"):
                continue

            summary.assimilate(Row(
                chr=row["#chr"],
                start=int(row["start"]),
                end=int(row["end"]),
                name=row["name"],
                beta_est=row["beta_est"],
                strand=row["strand"],
                unmod=int(row["unmod"]),
                mod=int(row["mod"]),
                no_snp=int(row["no_snp"]),
                snp=int(row["snp"]),
                coverage=int(row["coverage"]),
                genotype=row["genotype"],
                gt_p_score=int(row["gt_p_score"]),
                gt_conf_score=int(row["gt_conf_score"]),
            ))

    with open(output_path, "w") as out_fobj:
        summary.write(out_fobj)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Summarize mod/unmod counts and compute overall methylation rate."
    )
    parser.add_argument(
        "-i", "--input-file",
        required=True,
        help="Path to the rastair input mods file"
    )
    parser.add_argument(
        "-o", "--output-dir",
        default=None,
        help="Optional directory to write the summary output into. "
             "Output filename will match the input filename with '.summary'."
    )

    args = parser.parse_args()

    base = os.path.basename(args.input_file)
    summary_name = f"{os.path.splitext(base)[0]}.summary"

    if args.output_dir:
        os.makedirs(args.output_dir, exist_ok=True)
        summary_out = os.path.join(args.output_dir, summary_name)
    else:
        summary_out = summary_name

    main(args.input_file, summary_out)
