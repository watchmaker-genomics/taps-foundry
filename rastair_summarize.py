#!/usr/bin/env python3

import io
import csv
import argparse
from typing import NamedTuple
from dataclasses import dataclass


class Row(NamedTuple):
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
    total_mod: int = 0
    total_unmod: int = 0
    covered_positions: int = 0

    def assimilate(self, row: Row):
        self.total_mod += row.mod
        self.total_unmod += row.unmod

        if (row.mod + row.unmod) > 0:
            self.covered_positions += 1

    def write(self, out_fobj: io.TextIOWrapper):
        writer = csv.DictWriter(out_fobj, fieldnames=[
            "total_mod",
            "total_unmod",
            "methylation_rate",
            "covered_positions",
        ], delimiter="\t")
        methylation_rate = self.total_mod / (self.total_mod + self.total_unmod)
        writer.writeheader()
        writer.writerow({
            "total_mod": self.total_mod,
            "total_unmod": self.total_unmod,
            "methylation_rate": round(methylation_rate, 5),
            "covered_positions": self.covered_positions,
        })


def main(rastair_in: str, summary_out: str):
    summary = Summary()

    with open(rastair_in, newline="") as csvfile:
        reader = csv.DictReader(csvfile, delimiter="\t")
        for row in reader:
            row = Row(
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
            )
            summary.assimilate(row)

    with open(summary_out, "w") as out_fobj:
        summary.write(out_fobj)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--rastair_in', help='rastair result file to summarize')
    parser.add_argument('--summary_out', help='output file name for the summarized results')
    args = parser.parse_args()

    main(args.rastair_in, args.summary_out)
