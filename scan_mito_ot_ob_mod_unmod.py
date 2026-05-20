#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import pysam


@dataclass(frozen=True)
class TrimMask:
    r1_start: int
    r1_end: int
    r2_start: int
    r2_end: int

    @classmethod
    def parse(cls, value: str) -> "TrimMask":
        parts = [int(x) for x in value.split(",")]
        if len(parts) != 4:
            raise ValueError(f"Expected four comma-separated integers, got: {value}")
        return cls(*parts)

    def start_end_for_read(self, is_read1: bool) -> tuple[int, int]:
        if is_read1:
            return self.r1_start, self.r1_end
        return self.r2_start, self.r2_end


@dataclass
class Summary:
    group: str
    ref_base: str
    mod_base: str
    unmod_base: str
    positions: int = 0
    covered_positions: int = 0
    mod: int = 0
    unmod: int = 0
    other: int = 0

    @property
    def informative(self) -> int:
        return self.mod + self.unmod

    @property
    def depth(self) -> int:
        return self.mod + self.unmod + self.other

    @property
    def mod_fraction(self) -> float:
        if not self.informative:
            return 0.0
        return self.mod / self.informative


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Scan a mitochondrial contig for OT C->T and OB G->A observations, "
            "counting MOD and UNMOD base calls from a BAM."
        )
    )
    parser.add_argument("--bam", type=Path, required=True, help="Sorted and indexed BAM.")
    parser.add_argument("--fasta", type=Path, required=True, help="Indexed reference FASTA.")
    parser.add_argument("--contig", default="chrM", help="Contig to scan.")
    parser.add_argument("--min-mapq", type=int, default=20, help="Minimum mapping quality.")
    parser.add_argument("--min-baseq", type=int, default=30, help="Minimum base quality.")
    parser.add_argument(
        "--include-flags",
        type=int,
        default=3,
        help="Require all these SAM flag bits, matching rastair defaults.",
    )
    parser.add_argument(
        "--exclude-flags",
        type=int,
        default=3852,
        help="Exclude reads with any of these SAM flag bits, matching rastair defaults.",
    )
    parser.add_argument(
        "--nOT",
        default="0,0,0,0",
        help="OT trim mask as r1_start,r1_end,r2_start,r2_end.",
    )
    parser.add_argument(
        "--nOB",
        default="0,0,0,0",
        help="OB trim mask as r1_start,r1_end,r2_start,r2_end.",
    )
    parser.add_argument(
        "--max-depth",
        type=int,
        default=200000,
        help="Maximum pileup depth per position.",
    )
    parser.add_argument(
        "--out-tsv",
        type=Path,
        required=True,
        help="Per-position output TSV.",
    )
    parser.add_argument(
        "--summary-tsv",
        type=Path,
        required=True,
        help="Summary TSV.",
    )
    return parser.parse_args()


def classify_group(read: pysam.AlignedSegment) -> str | None:
    # Directional paired-end convention used by rastair's OT/OB masks:
    # OT: R1 forward, R2 reverse
    # OB: R1 reverse, R2 forward
    if read.is_read1 and not read.is_reverse:
        return "OT"
    if read.is_read2 and read.is_reverse:
        return "OT"
    if read.is_read1 and read.is_reverse:
        return "OB"
    if read.is_read2 and not read.is_reverse:
        return "OB"
    return None


def read_passes(read: pysam.AlignedSegment, include_flags: int, exclude_flags: int, min_mapq: int) -> bool:
    if (read.flag & include_flags) != include_flags:
        return False
    if read.flag & exclude_flags:
        return False
    if read.mapping_quality < min_mapq:
        return False
    return True


def base_is_trimmed(
    read: pysam.AlignedSegment,
    query_pos: int,
    group: str,
    nOT: TrimMask,
    nOB: TrimMask,
) -> bool:
    query_len = read.query_length
    if query_len is None:
        return True
    mask = nOT if group == "OT" else nOB
    start_trim, end_trim = mask.start_end_for_read(read.is_read1)
    if query_pos < start_trim:
        return True
    if end_trim and query_pos >= query_len - end_trim:
        return True
    return False


def write_summary(path: Path, rows: list[Summary]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "group",
                "ref_base",
                "mod_base",
                "unmod_base",
                "positions",
                "covered_positions",
                "mod",
                "unmod",
                "other",
                "informative",
                "depth",
                "mod_fraction",
            ]
        )
        for row in rows:
            writer.writerow(
                [
                    row.group,
                    row.ref_base,
                    row.mod_base,
                    row.unmod_base,
                    row.positions,
                    row.covered_positions,
                    row.mod,
                    row.unmod,
                    row.other,
                    row.informative,
                    row.depth,
                    f"{row.mod_fraction:.12f}",
                ]
            )


def main() -> None:
    args = parse_args()
    nOT = TrimMask.parse(args.nOT)
    nOB = TrimMask.parse(args.nOB)
    args.out_tsv.parent.mkdir(parents=True, exist_ok=True)
    args.summary_tsv.parent.mkdir(parents=True, exist_ok=True)

    summary_by_group = {
        "OT": Summary(group="OT", ref_base="C", mod_base="T", unmod_base="C"),
        "OB": Summary(group="OB", ref_base="G", mod_base="A", unmod_base="G"),
    }

    with pysam.AlignmentFile(args.bam, "rb") as bam, pysam.FastaFile(str(args.fasta)) as fasta, args.out_tsv.open(
        "w", newline=""
    ) as out_handle:
        contig_len = fasta.get_reference_length(args.contig)
        writer = csv.writer(out_handle, delimiter="\t")
        writer.writerow(
            [
                "chrom",
                "pos_1based",
                "group",
                "ref_base",
                "mod_base",
                "unmod_base",
                "mod",
                "unmod",
                "other",
                "informative",
                "depth",
                "mod_fraction",
            ]
        )

        for pileup_column in bam.pileup(
            contig=args.contig,
            start=0,
            stop=contig_len,
            truncate=True,
            stepper="all",
            fastafile=fasta,
            max_depth=args.max_depth,
            min_base_quality=0,
            ignore_overlaps=False,
        ):
            pos0 = pileup_column.reference_pos
            ref_base = fasta.fetch(args.contig, pos0, pos0 + 1).upper()
            if ref_base == "C":
                group = "OT"
                mod_base = "T"
                unmod_base = "C"
            elif ref_base == "G":
                group = "OB"
                mod_base = "A"
                unmod_base = "G"
            else:
                continue

            summary = summary_by_group[group]
            summary.positions += 1

            mod = 0
            unmod = 0
            other = 0

            for pileup_read in pileup_column.pileups:
                if pileup_read.is_del or pileup_read.is_refskip:
                    continue
                read = pileup_read.alignment
                if not read_passes(read, args.include_flags, args.exclude_flags, args.min_mapq):
                    continue

                read_group = classify_group(read)
                if read_group != group:
                    continue

                query_pos = pileup_read.query_position
                if query_pos is None:
                    continue
                if base_is_trimmed(read, query_pos, group, nOT, nOB):
                    continue

                qualities = read.query_qualities
                if qualities is None or qualities[query_pos] < args.min_baseq:
                    continue

                base = read.query_sequence[query_pos].upper()
                if base == mod_base:
                    mod += 1
                elif base == unmod_base:
                    unmod += 1
                else:
                    other += 1

            informative = mod + unmod
            depth = informative + other
            if depth:
                summary.covered_positions += 1
                summary.mod += mod
                summary.unmod += unmod
                summary.other += other
                mod_fraction = mod / informative if informative else 0.0
                writer.writerow(
                    [
                        args.contig,
                        pos0 + 1,
                        group,
                        ref_base,
                        mod_base,
                        unmod_base,
                        mod,
                        unmod,
                        other,
                        informative,
                        depth,
                        f"{mod_fraction:.12f}",
                    ]
                )

    combined = Summary(group="combined", ref_base="C/G", mod_base="T/A", unmod_base="C/G")
    for row in summary_by_group.values():
        combined.positions += row.positions
        combined.covered_positions += row.covered_positions
        combined.mod += row.mod
        combined.unmod += row.unmod
        combined.other += row.other

    write_summary(args.summary_tsv, [summary_by_group["OT"], summary_by_group["OB"], combined])


if __name__ == "__main__":
    main()
