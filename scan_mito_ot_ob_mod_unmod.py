#!/usr/bin/env python3

"""
Scan a mitochondrial contig for OT C->T and OB G->A base-call observations
and count MOD vs UNMOD reads at each position from an aligned BAM.
"""

from __future__ import annotations

import argparse
import csv
import sys
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import NoReturn

try:
    import pysam
except ImportError:
    sys.stderr.write(
        "Error: pysam is required but not installed. Install with `pip install pysam`.\n"
    )
    sys.exit(1)


def die(message: str) -> NoReturn:
    """
    Print an error message to stderr and exit with status 1.

    @param message: Human-readable error description (no trailing newline).
    """
    sys.stderr.write(f"Error: {message}\n")
    sys.exit(1)


class Group(Enum):
    """
    Read-orientation group used by TAPS/rastair-style methylation calling.

    OT = Original Top strand: the stranded read is forward (and its mate, if
        present, is reverse).
    OB = Original Bottom strand: the stranded read is reverse (and its mate,
        if present, is forward).
    """

    OT = "OT"
    OB = "OB"

    def __str__(self) -> str:
        return self.value


class StrandedRead(Enum):
    """
    Which mate of a paired-end library is the stranded (directional) read.

    R1 = SAM read 1 carries the original strand (standard directional libraries).
    R2 = SAM read 2 carries the original strand (e.g. PBAT / swapped libraries).

    Single-end reads ignore this setting: the single read is always the
    stranded read.
    """

    R1 = "R1"
    R2 = "R2"

    def __str__(self) -> str:
        return self.value


# ---------------------------------------------------------------------------
# Defaults — single source of truth for every CLI default and the TSV schema.
# Change here, the rest of the module follows.
# ---------------------------------------------------------------------------
DEFAULT_CONTIG: str = "chrM"
DEFAULT_MIN_MAPQ: int = 20
DEFAULT_MIN_BASEQ: int = 30
DEFAULT_INCLUDE_FLAGS: int = 3       # 0x1 paired + 0x2 proper-pair
DEFAULT_EXCLUDE_FLAGS: int = 3852    # unmapped/mate-unmapped/secondary/QC-fail/dup/supplementary
DEFAULT_TRIM_MASK: str = "0,0,0,0"   # r1_5p,r1_3p,r2_5p,r2_3p
DEFAULT_STRANDED_READ: StrandedRead = StrandedRead.R1
DEFAULT_MAX_DEPTH: int = 200_000
DEFAULT_MAX_SOFTCLIP: int = 5        # reject reads with >5 total soft-clipped bases

POSITION_TSV_HEADER: list[str] = [
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
    "nonmeth_mod",
    "nonmeth_unmod",
    "nonmeth_other",
    "nonmeth_informative",
    "nonmeth_depth",
    "nonmeth_mod_fraction",
    "beta_mod_fraction",
]


@dataclass(frozen=True)
class TrimMask:
    """
    Per-read trim mask in rastair v1 convention.

    @param r1_start: Bases to mask from the 5' end of R1 (or single-end).
    @param r1_end: Bases to mask from the 3' end of R1 (or single-end).
    @param r2_start: Bases to mask from the 5' end of R2.
    @param r2_end: Bases to mask from the 3' end of R2.
    """

    r1_start: int
    r1_end: int
    r2_start: int
    r2_end: int

    @classmethod
    def parse(cls, value: str) -> "TrimMask":
        """
        Parse a comma-separated trim mask string into a TrimMask.

        @param value: Four comma-separated integers as 'r1_start,r1_end,r2_start,r2_end'.
        @return: TrimMask populated from the parsed values.
        """
        parts = [int(x) for x in value.split(",")]
        if len(parts) != 4:
            raise ValueError(f"Expected four comma-separated integers, got: {value}")
        return cls(*parts)

    def start_end_for_read(self, is_read1: bool) -> tuple[int, int]:
        """
        Return the (5', 3') trim pair for a single read.

        @param is_read1: True for R1 (or single-end), False for R2.
        @return: Tuple of (start_trim, end_trim) where start is the 5' and end the 3' in the read's original orientation.
        """
        if is_read1:
            return self.r1_start, self.r1_end
        return self.r2_start, self.r2_end


@dataclass
class Summary:
    """
    Per-group aggregate counters across the scanned contig.

    @param group: Group label written to the summary TSV (e.g. 'OT', 'OB', 'combined').
    @param ref_base: Reference base for this group.
    @param mod_base: Read base interpreted as modified (e.g. 'T' for OT).
    @param unmod_base: Read base interpreted as unmodified (e.g. 'C' for OT).
    @param positions: Reference positions in this group that were scanned.
    @param covered_positions: Positions with at least one passing read.
    @param mod: Total modified base calls observed on the methylated strand.
    @param unmod: Total unmodified base calls observed on the methylated strand.
    @param other: Total base calls that were neither mod nor unmod on the methylated strand.
    @param nonmeth_mod: Same as mod, but tallied from opposite-strand reads at the position (control / background channel).
    @param nonmeth_unmod: Same as unmod, but tallied from opposite-strand reads.
    @param nonmeth_other: Same as other, but tallied from opposite-strand reads.
    """

    group: str
    ref_base: str
    mod_base: str
    unmod_base: str
    positions: int = 0
    covered_positions: int = 0
    mod: int = 0
    unmod: int = 0
    other: int = 0
    nonmeth_mod: int = 0
    nonmeth_unmod: int = 0
    nonmeth_other: int = 0

    @property
    def informative(self) -> int:
        """
        Mod + unmod calls — the denominator used for mod_fraction.

        @return: Total informative base calls.
        """
        return self.mod + self.unmod

    @property
    def depth(self) -> int:
        """
        Mod + unmod + other calls — total base calls accepted at this group.

        @return: Total base calls.
        """
        return self.mod + self.unmod + self.other

    @property
    def mod_fraction(self) -> float:
        """
        Fraction of informative calls that were modified.

        @return: mod / informative, or 0.0 when no informative calls exist.
        """
        if not self.informative:
            return 0.0
        return self.mod / self.informative

    @property
    def nonmeth_informative(self) -> int:
        """
        Non-meth mod + unmod calls — denominator for nonmeth_mod_fraction.

        @return: Total informative non-meth base calls.
        """
        return self.nonmeth_mod + self.nonmeth_unmod

    @property
    def nonmeth_depth(self) -> int:
        """
        Non-meth mod + unmod + other calls — total non-meth base calls.

        @return: Total non-meth base calls.
        """
        return self.nonmeth_mod + self.nonmeth_unmod + self.nonmeth_other

    @property
    def nonmeth_mod_fraction(self) -> float:
        """
        Fraction of informative non-meth calls that match the mod base.

        Expected to be ~sequencing-error rate; meaningfully elevated values
        flag candidate sites of real variation rather than TAPS signal.

        @return: nonmeth_mod / nonmeth_informative, or 0.0 when empty.
        """
        if not self.nonmeth_informative:
            return 0.0
        return self.nonmeth_mod / self.nonmeth_informative

    @property
    def beta_mod_fraction(self) -> float:
        """
        Error-corrected methylation estimate using the non-meth strand as the
        sequencing/deamination error rate:

            beta = (mod_fraction - nonmeth_mod_fraction) / (1 - nonmeth_mod_fraction)

        Clamped to >= 0. Falls back to mod_fraction when no non-meth reads
        are available (no correction possible).

        @return: Error-corrected mod fraction in [0, 1].
        """
        obs = self.mod_fraction
        if not self.nonmeth_informative:
            return obs
        err = self.nonmeth_mod_fraction
        if err >= 1.0:
            return obs
        corrected = (obs - err) / (1.0 - err)
        return max(0.0, corrected)


@dataclass
class PositionCounts:
    """
    Per-position tally of meth- and non-meth-strand base calls after dedup.

    Mirrors the count math on `Summary` so per-position rows and aggregate
    rows expose the same derived metrics.

    @param mod: Meth-strand calls matching mod_base.
    @param unmod: Meth-strand calls matching unmod_base.
    @param other: Meth-strand calls matching neither.
    @param nonmeth_mod: Opposite-strand calls matching mod_base.
    @param nonmeth_unmod: Opposite-strand calls matching unmod_base.
    @param nonmeth_other: Opposite-strand calls matching neither.
    """

    mod: int = 0
    unmod: int = 0
    other: int = 0
    nonmeth_mod: int = 0
    nonmeth_unmod: int = 0
    nonmeth_other: int = 0

    @property
    def informative(self) -> int:
        return self.mod + self.unmod

    @property
    def depth(self) -> int:
        return self.informative + self.other

    @property
    def mod_fraction(self) -> float:
        return self.mod / self.informative if self.informative else 0.0

    @property
    def nonmeth_informative(self) -> int:
        return self.nonmeth_mod + self.nonmeth_unmod

    @property
    def nonmeth_depth(self) -> int:
        return self.nonmeth_informative + self.nonmeth_other

    @property
    def nonmeth_mod_fraction(self) -> float:
        return (
            self.nonmeth_mod / self.nonmeth_informative
            if self.nonmeth_informative
            else 0.0
        )

    @property
    def beta_mod_fraction(self) -> float:
        obs = self.mod_fraction
        if not self.nonmeth_informative:
            return obs
        err = self.nonmeth_mod_fraction
        if err >= 1.0:
            return obs
        return max(0.0, (obs - err) / (1.0 - err))


def parse_args() -> argparse.Namespace:
    """
    Parse command-line arguments for the scanner.

    @return: argparse.Namespace populated from sys.argv.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Scan a mitochondrial contig for OT C->T and OB G->A observations, "
            "counting MOD and UNMOD base calls from a BAM."
        )
    )
    parser.add_argument("--bam", type=Path, required=True, help="Sorted and indexed BAM.")
    parser.add_argument("--fasta", type=Path, required=True, help="Indexed reference FASTA.")
    parser.add_argument("--contig", default=DEFAULT_CONTIG, help=f"Contig to scan. Default: {DEFAULT_CONTIG}.")
    parser.add_argument("--min-mapq", type=int, default=DEFAULT_MIN_MAPQ, help=f"Minimum mapping quality. Default: {DEFAULT_MIN_MAPQ}.")
    parser.add_argument("--min-baseq", type=int, default=DEFAULT_MIN_BASEQ, help=f"Minimum base quality. Default: {DEFAULT_MIN_BASEQ}.")
    parser.add_argument(
        "--include-flags",
        type=int,
        default=DEFAULT_INCLUDE_FLAGS,
        help=(
            f"Require all these SAM flag bits (rastair default: {DEFAULT_INCLUDE_FLAGS} = 0x1 paired + 0x2 proper-pair). "
            "Use 0 for single-end data."
        ),
    )
    parser.add_argument(
        "--exclude-flags",
        type=int,
        default=DEFAULT_EXCLUDE_FLAGS,
        help=(
            f"Exclude reads with any of these SAM flag bits (rastair default: {DEFAULT_EXCLUDE_FLAGS} = "
            "0x4 unmapped + 0x8 mate-unmapped + 0x100 secondary + 0x200 vendor-QC-fail + "
            "0x400 PCR/optical duplicate + 0x800 supplementary). "
            "Drops unmapped/orphan reads, non-primary alignments, QC failures, and dupes."
        ),
    )
    parser.add_argument(
        "--nOT",
        default=DEFAULT_TRIM_MASK,
        help="OT trim mask as r1_start,r1_end,r2_start,r2_end (r1/r2 = SAM read number, not the stranded-read setting).",
    )
    parser.add_argument(
        "--nOB",
        default=DEFAULT_TRIM_MASK,
        help="OB trim mask as r1_start,r1_end,r2_start,r2_end (r1/r2 = SAM read number, not the stranded-read setting).",
    )
    parser.add_argument(
        "--stranded-read",
        choices=[member.value for member in StrandedRead],
        default=DEFAULT_STRANDED_READ.value,
        help=f"Which paired-end mate carries the original DNA strand. Default: {DEFAULT_STRANDED_READ.value}. Ignored for single-end data.",
    )
    parser.add_argument(
        "--max-depth",
        type=int,
        default=DEFAULT_MAX_DEPTH,
        help=f"Maximum pileup depth per position. Default: {DEFAULT_MAX_DEPTH}.",
    )
    parser.add_argument(
        "--max-softclip",
        type=int,
        default=DEFAULT_MAX_SOFTCLIP,
        help=(
            f"Maximum total soft-clipped bases tolerated per read. Reads with more soft clip "
            f"are dropped. Default: {DEFAULT_MAX_SOFTCLIP} (i.e. reject reads with {DEFAULT_MAX_SOFTCLIP + 1} or more soft-clipped bases). "
            "Heavy soft-clipping is associated with adapter contamination and "
            "mis-alignment at fragment ends, both of which can bias methylation calls."
        ),
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


def classify_group(read: pysam.AlignedSegment, stranded_read: StrandedRead) -> Group:
    """
    Decide whether a read belongs to the OT or OB strand group.

    OT = stranded read forward (mate, if any, is reverse).
    OB = stranded read reverse (mate, if any, is forward).

    @param read: Aligned read to classify.
    @param stranded_read: Which mate of a paired-end library is the stranded read. Single-end reads ignore this and are always treated as the stranded read.
    @return: Group.OT or Group.OB.
    """
    if not read.is_paired:
        is_stranded = True
    elif stranded_read is StrandedRead.R1:
        is_stranded = read.is_read1
    else:
        is_stranded = read.is_read2

    forward = not read.is_reverse
    if is_stranded == forward:
        return Group.OT
    return Group.OB


def soft_clip_total(read: pysam.AlignedSegment) -> int:
    """
    Total bases soft-clipped across both ends of the read.

    Computed as `query_length - query_alignment_length`. `query_length`
    includes soft-clipped bases; `query_alignment_length` does not.

    @param read: Aligned read.
    @return: Total soft-clipped base count, or 0 when length info is missing.
    """
    qlen = read.query_length
    aqlen = read.query_alignment_length
    if qlen is None or aqlen is None:
        return 0
    return qlen - aqlen


def read_passes(
    read: pysam.AlignedSegment,
    include_flags: int,
    exclude_flags: int,
    max_softclip: int,
) -> bool:
    """
    Check whether a read passes the SAM flag and soft-clip filters.

    Mapping-quality filtering is handled by `pileup(min_mapping_quality=...)`
    upstream so it never reaches this function.

    @param read: Aligned read to test.
    @param include_flags: All bits required to be set in the read's SAM flag.
    @param exclude_flags: Bits that must not be set in the read's SAM flag.
    @param max_softclip: Maximum total soft-clipped bases the read may have.
        Reads with more soft clip than this are rejected.
    @return: True if the read should be kept, False otherwise.
    """
    if (read.flag & include_flags) != include_flags:
        return False
    if read.flag & exclude_flags:
        return False
    if soft_clip_total(read) > max_softclip:
        return False
    return True


def base_is_trimmed(
    read: pysam.AlignedSegment,
    query_pos: int,
    group: Group,
    nOT: TrimMask,
    nOB: TrimMask,
) -> bool:
    """
    Check whether a base call falls inside the masked 5'/3' region of its read.

    Trim values are applied in read orientation (rastair v1 convention), so for
    reverse-strand reads the start/end indices are swapped before being compared
    against query_position.

    @param read: Aligned read containing the base call.
    @param query_pos: Position of the base call in read.query_sequence.
    @param group: OT or OB classification of the read.
    @param nOT: Trim mask for OT reads.
    @param nOB: Trim mask for OB reads.
    @return: True if the base falls inside the trim region and should be skipped.
    """
    query_len = read.query_length
    if query_len is None:
        return True
    mask = nOT if group is Group.OT else nOB
    start_trim, end_trim = mask.start_end_for_read(not read.is_read2)
    # `start` = 5' of the original read, `end` = 3' (rastair v1 convention).
    # BAM stores SEQ reference-oriented, so for reverse-strand reads the 5'
    # end is at the right of query_sequence — swap the two trim values.
    if read.is_reverse:
        start_trim, end_trim = end_trim, start_trim
    if query_pos < start_trim:
        return True
    if end_trim and query_pos >= query_len - end_trim:
        return True
    return False


def write_summary(path: Path, rows: list[Summary]) -> None:
    """
    Write the per-group summary TSV.

    @param path: Destination TSV path. Parent directory must already exist.
    @param rows: Summary rows to write, one per output line.
    """
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
                "nonmeth_mod",
                "nonmeth_unmod",
                "nonmeth_other",
                "nonmeth_informative",
                "nonmeth_depth",
                "nonmeth_mod_fraction",
                "beta_mod_fraction",
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
                    row.nonmeth_mod,
                    row.nonmeth_unmod,
                    row.nonmeth_other,
                    row.nonmeth_informative,
                    row.nonmeth_depth,
                    f"{row.nonmeth_mod_fraction:.12f}",
                    f"{row.beta_mod_fraction:.12f}",
                ]
            )


def validate_input_paths(args: argparse.Namespace) -> None:
    """
    Check that the BAM, FASTA, and FASTA index exist on disk.

    @param args: Parsed CLI namespace with `bam` and `fasta` Path fields.
    """
    if not args.bam.is_file():
        die(f"BAM file not found: {args.bam}")
    if not args.fasta.is_file():
        die(f"FASTA file not found: {args.fasta}")
    fasta_index = Path(f"{args.fasta}.fai")
    if not fasta_index.is_file():
        die(
            f"FASTA index not found: {fasta_index}. "
            f"Generate it with `samtools faidx {args.fasta}`."
        )


def validate_bam_and_contig(
    bam: pysam.AlignmentFile, fasta: pysam.FastaFile, args: argparse.Namespace
) -> None:
    """
    Check that the BAM is indexed and that the chosen contig exists in both
    the BAM and the FASTA.

    @param bam: Opened BAM handle.
    @param fasta: Opened FASTA handle.
    @param args: Parsed CLI namespace; uses `bam`, `fasta`, and `contig`.
    """
    if not bam.has_index():
        die(
            f"BAM file is not indexed: {args.bam}. "
            f"Generate an index with `samtools index {args.bam}`."
        )
    for name, refs, src in (
        ("BAM", bam.references, args.bam),
        ("FASTA", fasta.references, args.fasta),
    ):
        if args.contig not in refs:
            preview = ", ".join(refs[:20])
            suffix = ", ..." if len(refs) > 20 else ""
            die(
                f"Contig '{args.contig}' is not present in {src}. "
                f"Available contigs: {preview}{suffix}"
            )


def count_position(
    pileup_column: "pysam.PileupColumn",
    position_group: Group,
    mod_base: str,
    unmod_base: str,
    args: argparse.Namespace,
    stranded_read: StrandedRead,
    nOT: TrimMask,
    nOB: TrimMask,
) -> PositionCounts:
    """
    Walk a pileup column, apply filters, dedup overlapping mates by query name
    (keeping the higher-BQ call), and return separate meth- and non-meth-strand
    tallies.

    The non-meth dict captures opposite-strand reads at the same position as
    a background-rate / variation channel.

    @param pileup_column: Pileup column from `bam.pileup`.
    @param position_group: OT for reference-C positions, OB for reference-G.
    @param mod_base: The base call interpreted as modified at this position.
    @param unmod_base: The base call interpreted as unmodified.
    @param args: Parsed CLI namespace; uses `include_flags`, `exclude_flags`, `max_softclip`.
    @param stranded_read: Which mate carries the original strand.
    @param nOT: Trim mask for OT reads.
    @param nOB: Trim mask for OB reads.
    @return: PositionCounts with meth- and non-meth-strand tallies.
    """
    candidates_meth: dict[str, tuple[str, int]] = {}
    candidates_nonmeth: dict[str, tuple[str, int]] = {}

    for pileup_read in pileup_column.pileups:
        if pileup_read.is_del or pileup_read.is_refskip:
            continue
        read = pileup_read.alignment
        if not read_passes(read, args.include_flags, args.exclude_flags, args.max_softclip):
            continue

        read_group = classify_group(read, stranded_read)
        target = candidates_meth if read_group == position_group else candidates_nonmeth

        query_pos = pileup_read.query_position
        if query_pos is None:
            continue
        # Pass the read's own group, not the position's, so opposite-strand
        # (non-meth control) reads get trimmed with their own mask. Matters
        # only when nOT != nOB; preserves correctness of nonmeth counts and
        # therefore beta_mod_fraction.
        if base_is_trimmed(read, query_pos, read_group, nOT, nOB):
            continue

        # Defensive: pysam returns None when the BAM record's QUAL field is
        # `*` (no quality scores stored). htslib's min_base_quality on the
        # pileup engine already filters those reads when --min-baseq > 0,
        # but we double-check here so qualities[query_pos] can never KeyError.
        qualities = read.query_qualities
        if qualities is None:
            continue

        bq = qualities[query_pos]
        base = read.query_sequence[query_pos].upper()
        existing = target.get(read.query_name)
        if existing is None or bq > existing[1]:
            target[read.query_name] = (base, bq)

    counts = PositionCounts()
    for base, _bq in candidates_meth.values():
        if base == mod_base:
            counts.mod += 1
        elif base == unmod_base:
            counts.unmod += 1
        else:
            counts.other += 1
    for base, _bq in candidates_nonmeth.values():
        if base == mod_base:
            counts.nonmeth_mod += 1
        elif base == unmod_base:
            counts.nonmeth_unmod += 1
        else:
            counts.nonmeth_other += 1
    return counts


def write_position_row(
    writer: "csv._writer",
    contig: str,
    pos_1based: int,
    position_group: Group,
    ref_base: str,
    mod_base: str,
    unmod_base: str,
    counts: PositionCounts,
) -> None:
    """
    Append one row to the per-position TSV.

    @param writer: csv.writer for the per-position TSV.
    @param contig: Reference contig name.
    @param pos_1based: 1-based reference position.
    @param position_group: OT or OB classification of this position.
    @param ref_base: Reference base at the position.
    @param mod_base: Base call interpreted as modified.
    @param unmod_base: Base call interpreted as unmodified.
    @param counts: PositionCounts for this row.
    """
    writer.writerow(
        [
            contig,
            pos_1based,
            position_group,
            ref_base,
            mod_base,
            unmod_base,
            counts.mod,
            counts.unmod,
            counts.other,
            counts.informative,
            counts.depth,
            f"{counts.mod_fraction:.12f}",
            counts.nonmeth_mod,
            counts.nonmeth_unmod,
            counts.nonmeth_other,
            counts.nonmeth_informative,
            counts.nonmeth_depth,
            f"{counts.nonmeth_mod_fraction:.12f}",
            f"{counts.beta_mod_fraction:.12f}",
        ]
    )


def accumulate_position(summary: Summary, counts: PositionCounts) -> None:
    """
    Fold a position's counts into its group summary (covered, mod, unmod, etc.).

    Called only for positions that had at least one meth-strand read.

    @param summary: Group summary to update in place.
    @param counts: PositionCounts produced for this position.
    """
    summary.covered_positions += 1
    summary.mod += counts.mod
    summary.unmod += counts.unmod
    summary.other += counts.other
    summary.nonmeth_mod += counts.nonmeth_mod
    summary.nonmeth_unmod += counts.nonmeth_unmod
    summary.nonmeth_other += counts.nonmeth_other


def build_combined_summary(summary_by_group: dict[Group, Summary]) -> Summary:
    """
    Build the 'combined' summary row by summing OT and OB counters.

    @param summary_by_group: Mapping of Group → Summary.
    @return: A new Summary labelled 'combined'.
    """
    combined = Summary(group="combined", ref_base="C/G", mod_base="T/A", unmod_base="C/G")
    for row in summary_by_group.values():
        combined.positions += row.positions
        combined.covered_positions += row.covered_positions
        combined.mod += row.mod
        combined.unmod += row.unmod
        combined.other += row.other
        combined.nonmeth_mod += row.nonmeth_mod
        combined.nonmeth_unmod += row.nonmeth_unmod
        combined.nonmeth_other += row.nonmeth_other
    return combined


def main() -> None:
    """
    Entry point: parse arguments, scan the contig, and write the per-position
    and summary TSVs.
    """
    args = parse_args()
    validate_input_paths(args)

    try:
        nOT = TrimMask.parse(args.nOT)
    except ValueError as exc:
        die(f"Invalid --nOT value '{args.nOT}': {exc}")
    try:
        nOB = TrimMask.parse(args.nOB)
    except ValueError as exc:
        die(f"Invalid --nOB value '{args.nOB}': {exc}")
    stranded_read = StrandedRead(args.stranded_read)
    args.out_tsv.parent.mkdir(parents=True, exist_ok=True)
    args.summary_tsv.parent.mkdir(parents=True, exist_ok=True)

    summary_by_group = {
        Group.OT: Summary(group=str(Group.OT), ref_base="C", mod_base="T", unmod_base="C"),
        Group.OB: Summary(group=str(Group.OB), ref_base="G", mod_base="A", unmod_base="G"),
    }

    with pysam.AlignmentFile(args.bam, "rb") as bam, pysam.FastaFile(str(args.fasta)) as fasta, args.out_tsv.open(
        "w", newline=""
    ) as out_handle:
        validate_bam_and_contig(bam, fasta, args)

        contig_len = fasta.get_reference_length(args.contig)
        writer = csv.writer(out_handle, delimiter="\t")
        writer.writerow(POSITION_TSV_HEADER)

        for pileup_column in bam.pileup(
            contig=args.contig,
            start=0,
            stop=contig_len,
            truncate=True,
            stepper="all",
            fastafile=fasta,
            max_depth=args.max_depth,
            min_base_quality=args.min_baseq,
            min_mapping_quality=args.min_mapq,
            ignore_overlaps=False,
        ):
            pos0 = pileup_column.reference_pos
            ref_base = fasta.fetch(args.contig, pos0, pos0 + 1).upper()
            if ref_base == "C":
                position_group, mod_base, unmod_base = Group.OT, "T", "C"
            elif ref_base == "G":
                position_group, mod_base, unmod_base = Group.OB, "A", "G"
            else:
                continue

            summary = summary_by_group[position_group]
            summary.positions += 1

            counts = count_position(
                pileup_column,
                position_group,
                mod_base,
                unmod_base,
                args,
                stranded_read,
                nOT,
                nOB,
            )

            if counts.depth:
                accumulate_position(summary, counts)
                write_position_row(
                    writer,
                    args.contig,
                    pos0 + 1,
                    position_group,
                    ref_base,
                    mod_base,
                    unmod_base,
                    counts,
                )

    combined = build_combined_summary(summary_by_group)
    write_summary(
        args.summary_tsv,
        [summary_by_group[Group.OT], summary_by_group[Group.OB], combined],
    )


if __name__ == "__main__":
    main()
