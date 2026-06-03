"""
Unit tests for scan_mito_ot_ob_mod_unmod.py — Easy + Medium tier.

Run with either:
    python -m unittest tests.test_scan_mito -v
    python -m unittest discover -s tests -v

No third-party deps beyond pysam (which the script already requires).
"""

from __future__ import annotations

import csv
import io
import sys
import unittest
from dataclasses import dataclass, field
from pathlib import Path
from types import SimpleNamespace

# Make the repo root importable so `import scan_mito_ot_ob_mod_unmod` works
# whether tests are run from the repo root or anywhere else.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import scan_mito_ot_ob_mod_unmod as smom  # noqa: E402


@dataclass
class FakeRead:
    """Stand-in for pysam.AlignedSegment with only the attrs the script reads."""

    is_paired: bool = True
    is_read1: bool = True
    is_read2: bool = False
    is_reverse: bool = False
    flag: int = 3  # paired + proper-pair
    mapping_quality: int = 60
    query_length: int = 100
    query_name: str = "frag1"
    query_sequence: str = "A" * 100
    query_qualities: list = field(default_factory=lambda: [40] * 100)


# ---------------------------------------------------------------------------
# TrimMask
# ---------------------------------------------------------------------------

class TestTrimMaskParse(unittest.TestCase):
    def test_happy_path(self):
        m = smom.TrimMask.parse("1,2,3,4")
        self.assertEqual((m.r1_start, m.r1_end, m.r2_start, m.r2_end), (1, 2, 3, 4))

    def test_all_zero(self):
        m = smom.TrimMask.parse("0,0,0,0")
        self.assertEqual((m.r1_start, m.r1_end, m.r2_start, m.r2_end), (0, 0, 0, 0))

    def test_wrong_field_count(self):
        with self.assertRaises(ValueError):
            smom.TrimMask.parse("1,2,3")

    def test_non_integer(self):
        with self.assertRaises(ValueError):
            smom.TrimMask.parse("1,2,3,foo")


class TestTrimMaskStartEndForRead(unittest.TestCase):
    def test_r1(self):
        m = smom.TrimMask(1, 2, 3, 4)
        self.assertEqual(m.start_end_for_read(True), (1, 2))

    def test_r2(self):
        m = smom.TrimMask(1, 2, 3, 4)
        self.assertEqual(m.start_end_for_read(False), (3, 4))


# ---------------------------------------------------------------------------
# Enums
# ---------------------------------------------------------------------------

class TestEnums(unittest.TestCase):
    def test_group_str(self):
        self.assertEqual(str(smom.Group.OT), "OT")
        self.assertEqual(str(smom.Group.OB), "OB")

    def test_stranded_read_str(self):
        self.assertEqual(str(smom.StrandedRead.R1), "R1")
        self.assertEqual(str(smom.StrandedRead.R2), "R2")


# ---------------------------------------------------------------------------
# Summary properties
# ---------------------------------------------------------------------------

class TestSummaryProperties(unittest.TestCase):
    def _summary(self, **kwargs):
        return smom.Summary("OT", "C", "T", "C", **kwargs)

    def test_zero_informative(self):
        self.assertEqual(self._summary().mod_fraction, 0.0)
        self.assertEqual(self._summary().informative, 0)
        self.assertEqual(self._summary().depth, 0)

    def test_mod_fraction(self):
        self.assertEqual(self._summary(mod=50, unmod=50).mod_fraction, 0.5)

    def test_informative_depth(self):
        s = self._summary(mod=50, unmod=50, other=10)
        self.assertEqual(s.informative, 100)
        self.assertEqual(s.depth, 110)

    def test_nonmeth_mod_fraction(self):
        s = self._summary(nonmeth_mod=10, nonmeth_unmod=90)
        self.assertEqual(s.nonmeth_mod_fraction, 0.1)

    def test_nonmeth_informative_zero(self):
        self.assertEqual(self._summary().nonmeth_mod_fraction, 0.0)


class TestSummaryBeta(unittest.TestCase):
    """Edge cases for the beta correction."""

    def _summary(self, **kwargs):
        return smom.Summary("OT", "C", "T", "C", **kwargs)

    def test_no_nonmeth_falls_back_to_mod_fraction(self):
        s = self._summary(mod=50, unmod=50)
        self.assertEqual(s.beta_mod_fraction, 0.5)

    def test_basic_correction(self):
        # obs=0.5, err=0.1 -> (0.5-0.1)/0.9 = 0.4444...
        s = self._summary(mod=50, unmod=50, nonmeth_mod=10, nonmeth_unmod=90)
        self.assertAlmostEqual(s.beta_mod_fraction, 0.4 / 0.9, places=10)

    def test_noise_above_signal_clamps_to_zero(self):
        # obs=0.05, err=0.20 -> negative -> clamp 0
        s = self._summary(mod=5, unmod=95, nonmeth_mod=20, nonmeth_unmod=80)
        self.assertEqual(s.beta_mod_fraction, 0.0)

    def test_signal_equals_noise_is_zero(self):
        # obs=0.10, err=0.10 -> 0
        s = self._summary(mod=10, unmod=90, nonmeth_mod=10, nonmeth_unmod=90)
        self.assertAlmostEqual(s.beta_mod_fraction, 0.0, places=10)

    def test_homozygous_alt_guard(self):
        # err == 1.0 -> divide-by-zero guard fires, returns raw mod_fraction
        s = self._summary(mod=10, unmod=0, nonmeth_mod=10, nonmeth_unmod=0)
        self.assertEqual(s.beta_mod_fraction, 1.0)

    def test_zero_mod_with_nonmeth_clamps_to_zero(self):
        # obs=0, err=0.1 -> (0 - 0.1)/0.9 = negative -> clamp 0
        s = self._summary(mod=0, unmod=100, nonmeth_mod=10, nonmeth_unmod=90)
        self.assertEqual(s.beta_mod_fraction, 0.0)


# ---------------------------------------------------------------------------
# PositionCounts properties — mirrors Summary math
# ---------------------------------------------------------------------------

class TestPositionCounts(unittest.TestCase):
    def test_mod_fraction(self):
        c = smom.PositionCounts(mod=3, unmod=7)
        self.assertAlmostEqual(c.mod_fraction, 0.3, places=10)

    def test_zero_informative(self):
        c = smom.PositionCounts()
        self.assertEqual(c.mod_fraction, 0.0)
        self.assertEqual(c.beta_mod_fraction, 0.0)

    def test_beta_full_correction(self):
        c = smom.PositionCounts(mod=50, unmod=50, nonmeth_mod=10, nonmeth_unmod=90)
        self.assertAlmostEqual(c.beta_mod_fraction, 0.4 / 0.9, places=10)

    def test_beta_homozygous_alt_guard(self):
        c = smom.PositionCounts(mod=10, unmod=0, nonmeth_mod=10, nonmeth_unmod=0)
        self.assertEqual(c.beta_mod_fraction, 1.0)

    def test_depth_includes_other(self):
        c = smom.PositionCounts(mod=1, unmod=1, other=3, nonmeth_other=5)
        self.assertEqual(c.depth, 5)
        self.assertEqual(c.nonmeth_depth, 5)


# ---------------------------------------------------------------------------
# classify_group
# ---------------------------------------------------------------------------

class TestClassifyGroupR1Stranded(unittest.TestCase):
    """stranded_read=R1 (default): R1 (or single-end) carries the strand."""

    def test_single_end_forward_is_ot(self):
        read = FakeRead(is_paired=False, is_read1=False, is_read2=False, is_reverse=False)
        self.assertEqual(smom.classify_group(read, smom.StrandedRead.R1), smom.Group.OT)

    def test_single_end_reverse_is_ob(self):
        read = FakeRead(is_paired=False, is_read1=False, is_read2=False, is_reverse=True)
        self.assertEqual(smom.classify_group(read, smom.StrandedRead.R1), smom.Group.OB)

    def test_paired_r1_forward_is_ot(self):
        read = FakeRead(is_paired=True, is_read1=True, is_read2=False, is_reverse=False)
        self.assertEqual(smom.classify_group(read, smom.StrandedRead.R1), smom.Group.OT)

    def test_paired_r2_reverse_is_ot(self):
        read = FakeRead(is_paired=True, is_read1=False, is_read2=True, is_reverse=True)
        self.assertEqual(smom.classify_group(read, smom.StrandedRead.R1), smom.Group.OT)

    def test_paired_r1_reverse_is_ob(self):
        read = FakeRead(is_paired=True, is_read1=True, is_read2=False, is_reverse=True)
        self.assertEqual(smom.classify_group(read, smom.StrandedRead.R1), smom.Group.OB)

    def test_paired_r2_forward_is_ob(self):
        read = FakeRead(is_paired=True, is_read1=False, is_read2=True, is_reverse=False)
        self.assertEqual(smom.classify_group(read, smom.StrandedRead.R1), smom.Group.OB)


class TestClassifyGroupR2Stranded(unittest.TestCase):
    """stranded_read=R2: R2 carries the strand (R1 is the mate)."""

    def test_paired_r2_forward_is_ot(self):
        read = FakeRead(is_paired=True, is_read1=False, is_read2=True, is_reverse=False)
        self.assertEqual(smom.classify_group(read, smom.StrandedRead.R2), smom.Group.OT)

    def test_paired_r1_forward_is_ob(self):
        read = FakeRead(is_paired=True, is_read1=True, is_read2=False, is_reverse=False)
        self.assertEqual(smom.classify_group(read, smom.StrandedRead.R2), smom.Group.OB)

    def test_single_end_unaffected(self):
        # Single-end ignores the stranded_read setting (no mate to flip with).
        read = FakeRead(is_paired=False, is_read1=False, is_read2=False, is_reverse=False)
        self.assertEqual(smom.classify_group(read, smom.StrandedRead.R2), smom.Group.OT)


# ---------------------------------------------------------------------------
# read_passes
# ---------------------------------------------------------------------------

class TestReadPasses(unittest.TestCase):
    def test_passes_defaults(self):
        read = FakeRead(flag=3)
        self.assertTrue(smom.read_passes(read, include_flags=3, exclude_flags=3852))

    def test_missing_required_bit(self):
        # flag=1 (paired only) — missing proper-pair bit
        read = FakeRead(flag=1)
        self.assertFalse(smom.read_passes(read, include_flags=3, exclude_flags=3852))

    def test_excluded_duplicate(self):
        # 0x400 = PCR/optical duplicate, in default exclude_flags=3852
        read = FakeRead(flag=3 | 0x400)
        self.assertFalse(smom.read_passes(read, include_flags=3, exclude_flags=3852))

    def test_single_end_with_include_zero(self):
        # Unpaired read passes when include_flags=0 (the recommended single-end setting)
        read = FakeRead(is_paired=False, flag=0)
        self.assertTrue(smom.read_passes(read, include_flags=0, exclude_flags=3852))


# ---------------------------------------------------------------------------
# base_is_trimmed — read-orientation semantics (rastair v1)
# ---------------------------------------------------------------------------

class TestBaseIsTrimmed(unittest.TestCase):
    def test_no_trim_returns_false(self):
        read = FakeRead(is_reverse=False, query_length=100)
        m = smom.TrimMask(0, 0, 0, 0)
        self.assertFalse(smom.base_is_trimmed(read, 50, smom.Group.OT, m, m))

    def test_forward_start_trim(self):
        read = FakeRead(is_reverse=False, query_length=100)
        m = smom.TrimMask(10, 0, 0, 0)  # r1_start = 10 (5' on a forward read)
        self.assertTrue(smom.base_is_trimmed(read, 5, smom.Group.OT, m, m))
        self.assertFalse(smom.base_is_trimmed(read, 10, smom.Group.OT, m, m))  # boundary

    def test_forward_end_trim(self):
        read = FakeRead(is_reverse=False, query_length=100)
        m = smom.TrimMask(0, 10, 0, 0)  # r1_end = 10 (3' on a forward read)
        self.assertTrue(smom.base_is_trimmed(read, 90, smom.Group.OT, m, m))
        self.assertFalse(smom.base_is_trimmed(read, 89, smom.Group.OT, m, m))  # boundary

    def test_reverse_swaps_5p_and_3p(self):
        """On a reverse-strand read the 5' end is at the right of query_sequence,
        so r1_start (5') applies near the right and r1_end (3') near the left."""
        read = FakeRead(is_reverse=True, query_length=100)
        m = smom.TrimMask(10, 0, 0, 0)  # 5'=10, 3'=0
        # After the strand swap, the 10-base trim falls on positions >= 90.
        self.assertTrue(smom.base_is_trimmed(read, 95, smom.Group.OT, m, m))
        self.assertFalse(smom.base_is_trimmed(read, 50, smom.Group.OT, m, m))
        # Left side has no 3' trim → not trimmed.
        self.assertFalse(smom.base_is_trimmed(read, 5, smom.Group.OT, m, m))

    def test_picks_r2_mask_for_r2_reads(self):
        # R2 read; the function uses r2_* fields.
        read = FakeRead(is_read1=False, is_read2=True, is_reverse=True, query_length=100)
        # r2_start=10 (5' of R2). R2 is reverse here, so 5' is at the right of SEQ.
        m = smom.TrimMask(0, 0, 10, 0)
        self.assertTrue(smom.base_is_trimmed(read, 95, smom.Group.OT, m, m))
        self.assertFalse(smom.base_is_trimmed(read, 50, smom.Group.OT, m, m))

    def test_missing_query_length_is_trimmed(self):
        read = FakeRead(query_length=None)
        m = smom.TrimMask(0, 0, 0, 0)
        self.assertTrue(smom.base_is_trimmed(read, 0, smom.Group.OT, m, m))


# ---------------------------------------------------------------------------
# accumulate_position
# ---------------------------------------------------------------------------

class TestAccumulatePosition(unittest.TestCase):
    def test_increments_covered_and_counts(self):
        s = smom.Summary("OT", "C", "T", "C")
        c = smom.PositionCounts(mod=3, unmod=7, other=1, nonmeth_mod=1, nonmeth_unmod=9, nonmeth_other=2)
        smom.accumulate_position(s, c)
        self.assertEqual(s.covered_positions, 1)
        self.assertEqual(s.mod, 3)
        self.assertEqual(s.unmod, 7)
        self.assertEqual(s.other, 1)
        self.assertEqual(s.nonmeth_mod, 1)
        self.assertEqual(s.nonmeth_unmod, 9)
        self.assertEqual(s.nonmeth_other, 2)

    def test_repeated_accumulation_sums(self):
        s = smom.Summary("OT", "C", "T", "C")
        c1 = smom.PositionCounts(mod=2, unmod=3)
        c2 = smom.PositionCounts(mod=4, unmod=1, nonmeth_mod=1)
        smom.accumulate_position(s, c1)
        smom.accumulate_position(s, c2)
        self.assertEqual(s.covered_positions, 2)
        self.assertEqual(s.mod, 6)
        self.assertEqual(s.unmod, 4)
        self.assertEqual(s.nonmeth_mod, 1)

    def test_does_not_touch_positions(self):
        # `positions` is incremented in main's pileup loop, not by accumulate_position.
        s = smom.Summary("OT", "C", "T", "C", positions=10)
        c = smom.PositionCounts(mod=1)
        smom.accumulate_position(s, c)
        self.assertEqual(s.positions, 10)
        self.assertEqual(s.covered_positions, 1)


# ---------------------------------------------------------------------------
# build_combined_summary
# ---------------------------------------------------------------------------

class TestBuildCombinedSummary(unittest.TestCase):
    def test_sums_ot_and_ob(self):
        ot = smom.Summary(
            "OT", "C", "T", "C",
            positions=10, covered_positions=8, mod=5, unmod=3, other=1,
            nonmeth_mod=1, nonmeth_unmod=2, nonmeth_other=3,
        )
        ob = smom.Summary(
            "OB", "G", "A", "G",
            positions=20, covered_positions=15, mod=10, unmod=4, other=2,
            nonmeth_mod=2, nonmeth_unmod=4, nonmeth_other=6,
        )
        combined = smom.build_combined_summary({smom.Group.OT: ot, smom.Group.OB: ob})
        self.assertEqual(combined.group, "combined")
        self.assertEqual(combined.ref_base, "C/G")
        self.assertEqual(combined.mod_base, "T/A")
        self.assertEqual(combined.positions, 30)
        self.assertEqual(combined.covered_positions, 23)
        self.assertEqual(combined.mod, 15)
        self.assertEqual(combined.unmod, 7)
        self.assertEqual(combined.other, 3)
        self.assertEqual(combined.nonmeth_mod, 3)
        self.assertEqual(combined.nonmeth_unmod, 6)
        self.assertEqual(combined.nonmeth_other, 9)


# ---------------------------------------------------------------------------
# write_position_row
# ---------------------------------------------------------------------------

class TestWritePositionRow(unittest.TestCase):
    def test_columns_and_values(self):
        buf = io.StringIO()
        writer = csv.writer(buf, delimiter="\t")
        counts = smom.PositionCounts(
            mod=3, unmod=7, other=0,
            nonmeth_mod=0, nonmeth_unmod=10, nonmeth_other=0,
        )
        smom.write_position_row(
            writer, "chrM", 100, smom.Group.OT, "C", "T", "C", counts
        )
        fields = buf.getvalue().rstrip("\r\n").split("\t")
        # 19 columns total
        self.assertEqual(len(fields), len(smom.POSITION_TSV_HEADER))
        # Spot-check a few high-value positions
        self.assertEqual(fields[0], "chrM")
        self.assertEqual(fields[1], "100")
        self.assertEqual(fields[2], "OT")
        self.assertEqual(fields[6], "3")  # mod
        self.assertEqual(fields[7], "7")  # unmod
        self.assertEqual(fields[11], f"{0.3:.12f}")  # mod_fraction
        # nonmeth: no mod calls so nonmeth_mod_fraction=0.0; beta falls back to mod_fraction
        self.assertEqual(fields[17], f"{0.0:.12f}")
        self.assertEqual(fields[18], f"{0.3:.12f}")


# ---------------------------------------------------------------------------
# write_summary
# ---------------------------------------------------------------------------

class TestWriteSummary(unittest.TestCase):
    def test_round_trip(self):
        import tempfile
        ot = smom.Summary("OT", "C", "T", "C", positions=2, covered_positions=2, mod=1, unmod=1)
        ob = smom.Summary("OB", "G", "A", "G", positions=3, covered_positions=3, mod=0, unmod=3)
        combined = smom.build_combined_summary({smom.Group.OT: ot, smom.Group.OB: ob})
        with tempfile.TemporaryDirectory() as td:
            out = Path(td) / "summary.tsv"
            smom.write_summary(out, [ot, ob, combined])
            lines = out.read_text().splitlines()
        self.assertEqual(len(lines), 4)  # header + 3 rows
        header = lines[0].split("\t")
        self.assertEqual(header[0], "group")
        self.assertEqual(header[-1], "beta_mod_fraction")
        self.assertTrue(lines[1].startswith("OT\tC\tT\tC"))
        self.assertTrue(lines[2].startswith("OB\tG\tA\tG"))
        self.assertTrue(lines[3].startswith("combined\tC/G\tT/A\tC/G"))
        # OT had no non-meth reads -> beta == mod_fraction == 0.5
        ot_fields = lines[1].split("\t")
        self.assertEqual(ot_fields[-1], f"{0.5:.12f}")


# ---------------------------------------------------------------------------
# validate_input_paths
# ---------------------------------------------------------------------------

class TestValidateInputPaths(unittest.TestCase):
    def _args(self, bam: Path, fasta: Path) -> SimpleNamespace:
        return SimpleNamespace(bam=bam, fasta=fasta)

    def test_missing_bam_exits(self):
        import tempfile
        with tempfile.TemporaryDirectory() as td:
            args = self._args(Path(td) / "nope.bam", Path(td) / "nope.fa")
            with self.assertRaises(SystemExit):
                smom.validate_input_paths(args)

    def test_missing_fasta_exits(self):
        import tempfile
        with tempfile.TemporaryDirectory() as td:
            bam = Path(td) / "ok.bam"
            bam.touch()
            args = self._args(bam, Path(td) / "nope.fa")
            with self.assertRaises(SystemExit):
                smom.validate_input_paths(args)

    def test_missing_fai_exits(self):
        import tempfile
        with tempfile.TemporaryDirectory() as td:
            bam = Path(td) / "ok.bam"
            bam.touch()
            fasta = Path(td) / "ok.fa"
            fasta.touch()
            args = self._args(bam, fasta)
            with self.assertRaises(SystemExit):
                smom.validate_input_paths(args)

    def test_all_present_no_exit(self):
        import tempfile
        with tempfile.TemporaryDirectory() as td:
            bam = Path(td) / "ok.bam"
            bam.touch()
            fasta = Path(td) / "ok.fa"
            fasta.touch()
            (Path(td) / "ok.fa.fai").touch()
            args = self._args(bam, fasta)
            # Should not raise
            smom.validate_input_paths(args)


# ---------------------------------------------------------------------------
# count_position — inner loop with FakePileupRead / FakePileupColumn
# ---------------------------------------------------------------------------

def _make_read(
    name: str,
    base: str,
    bq: int,
    *,
    is_read2: bool = False,
    is_reverse: bool = False,
    query_pos: int = 50,
    query_length: int = 100,
    flag: int = 3,
    mapq: int = 60,
    drop_qualities: bool = False,
) -> FakeRead:
    """Build a FakeRead with `base` placed at `query_pos` in the SEQ."""
    seq = ["N"] * query_length
    seq[query_pos] = base
    read = FakeRead(
        is_paired=True,
        is_read1=not is_read2,
        is_read2=is_read2,
        is_reverse=is_reverse,
        flag=flag,
        mapping_quality=mapq,
        query_length=query_length,
        query_name=name,
        query_sequence="".join(seq),
        query_qualities=[bq] * query_length,
    )
    if drop_qualities:
        read.query_qualities = None
    return read


def _pread(read: FakeRead, *, query_pos: int = 50, is_del: bool = False, is_refskip: bool = False) -> SimpleNamespace:
    return SimpleNamespace(
        alignment=read,
        query_position=query_pos,
        is_del=is_del,
        is_refskip=is_refskip,
    )


def _pcol(preads: list) -> SimpleNamespace:
    return SimpleNamespace(pileups=preads, reference_pos=99)


class TestCountPosition(unittest.TestCase):
    def _args(self, **kw):
        defaults = dict(include_flags=3, exclude_flags=3852)
        defaults.update(kw)
        return SimpleNamespace(**defaults)

    def _no_trim(self):
        return smom.TrimMask(0, 0, 0, 0)

    def _call(self, col, *, position_group=smom.Group.OT, mod_base="T", unmod_base="C",
              args=None, stranded=smom.StrandedRead.R1, nOT=None, nOB=None):
        return smom.count_position(
            col,
            position_group,
            mod_base,
            unmod_base,
            args or self._args(),
            stranded,
            nOT or self._no_trim(),
            nOB or self._no_trim(),
        )

    def test_single_ot_read_with_mod_base(self):
        col = _pcol([_pread(_make_read("frag1", "T", 40))])
        counts = self._call(col)
        self.assertEqual(counts.mod, 1)
        self.assertEqual(counts.unmod, 0)
        self.assertEqual(counts.other, 0)
        self.assertEqual(counts.nonmeth_mod, 0)

    def test_single_ot_read_with_unmod_base(self):
        col = _pcol([_pread(_make_read("frag1", "C", 40))])
        counts = self._call(col)
        self.assertEqual(counts.mod, 0)
        self.assertEqual(counts.unmod, 1)

    def test_other_base_goes_to_other(self):
        # G is neither mod (T) nor unmod (C)
        col = _pcol([_pread(_make_read("frag1", "G", 40))])
        counts = self._call(col)
        self.assertEqual(counts.other, 1)
        self.assertEqual(counts.mod, 0)
        self.assertEqual(counts.unmod, 0)

    def test_opposite_strand_goes_to_nonmeth(self):
        # At a C position (OT), an R2-forward read classifies as OB → non-meth bucket.
        read = _make_read("frag1", "T", 40, is_read2=True, is_reverse=False)
        counts = self._call(_pcol([_pread(read)]))
        self.assertEqual(counts.mod, 0)
        self.assertEqual(counts.nonmeth_mod, 1)

    def test_dedup_keeps_higher_bq(self):
        # Same fragment, both mates classify as OT (R1 fwd + R2 rev). Different bases.
        r1 = _make_read("frag1", "T", 30, is_read2=False, is_reverse=False)  # mod, lower BQ
        r2 = _make_read("frag1", "C", 40, is_read2=True, is_reverse=True)    # unmod, higher BQ
        counts = self._call(_pcol([_pread(r1), _pread(r2)]))
        # Higher-BQ (C) wins → unmod=1, mod=0
        self.assertEqual(counts.mod, 0)
        self.assertEqual(counts.unmod, 1)

    def test_dedup_tie_keeps_first(self):
        # Equal BQ — `bq > existing[1]` is strict, so first wins.
        r1 = _make_read("frag1", "T", 40, is_read2=False, is_reverse=False)
        r2 = _make_read("frag1", "C", 40, is_read2=True, is_reverse=True)
        counts = self._call(_pcol([_pread(r1), _pread(r2)]))
        self.assertEqual(counts.mod, 1)
        self.assertEqual(counts.unmod, 0)

    def test_different_fragments_both_counted(self):
        counts = self._call(_pcol([
            _pread(_make_read("fragA", "T", 40)),
            _pread(_make_read("fragB", "C", 40)),
        ]))
        self.assertEqual(counts.mod, 1)
        self.assertEqual(counts.unmod, 1)

    def test_skips_is_del(self):
        col = _pcol([_pread(_make_read("frag1", "T", 40), is_del=True)])
        self.assertEqual(self._call(col).depth, 0)

    def test_skips_is_refskip(self):
        col = _pcol([_pread(_make_read("frag1", "T", 40), is_refskip=True)])
        self.assertEqual(self._call(col).depth, 0)

    def test_skips_failed_read_passes(self):
        # 0x400 (PCR dup) is in exclude_flags=3852
        read = _make_read("frag1", "T", 40, flag=3 | 0x400)
        self.assertEqual(self._call(_pcol([_pread(read)])).depth, 0)

    def test_skips_none_query_position(self):
        # pysam returns query_position=None for soft-clipped / deletion-edge cases.
        col = _pcol([_pread(_make_read("frag1", "T", 40), query_pos=None)])
        self.assertEqual(self._call(col).depth, 0)

    def test_skips_trimmed_base(self):
        # 5' trim of 60 bases — position 50 is inside the trim region.
        read = _make_read("frag1", "T", 40, query_pos=50)
        col = _pcol([_pread(read, query_pos=50)])
        trim = smom.TrimMask(60, 0, 0, 0)
        counts = self._call(col, nOT=trim, nOB=smom.TrimMask(0, 0, 0, 0))
        self.assertEqual(counts.depth, 0)

    def test_opposite_strand_uses_its_own_trim_mask(self):
        # Regression: with asymmetric trim, an OB read at a reference-C (OT)
        # position must be trimmed by nOB, not nOT. Earlier code passed the
        # position's group to base_is_trimmed, biasing the nonmeth control
        # channel whenever nOT != nOB.
        # OB classification = R2 forward; place the base at query_pos=50.
        read = _make_read("frag1", "T", 40, is_read2=True, is_reverse=False, query_pos=50)
        col = _pcol([_pread(read, query_pos=50)])
        # nOT has no trim. nOB has 60-base 5' trim on R2 — under the bug, the
        # position group (OT) wins so nOT applies and the read is NOT trimmed.
        nOT = smom.TrimMask(0, 0, 0, 0)
        nOB = smom.TrimMask(0, 0, 60, 0)
        counts = self._call(col, nOT=nOT, nOB=nOB)
        # With the fix, the OB read is trimmed by nOB and is not counted as
        # non-meth at this position.
        self.assertEqual(counts.nonmeth_mod, 0)
        self.assertEqual(counts.nonmeth_unmod, 0)
        self.assertEqual(counts.nonmeth_other, 0)
        # And the meth channel was already empty (no OT reads).
        self.assertEqual(counts.depth, 0)

    def test_skips_when_qualities_is_none(self):
        read = _make_read("frag1", "T", 40, drop_qualities=True)
        self.assertEqual(self._call(_pcol([_pread(read)])).depth, 0)

    def test_mixed_meth_and_nonmeth_at_same_position(self):
        # OT (R1 fwd) and OB (R2 fwd) at a reference-C position
        ot = _make_read("fragA", "T", 40, is_read2=False, is_reverse=False)
        ob = _make_read("fragB", "C", 40, is_read2=True, is_reverse=False)
        counts = self._call(_pcol([_pread(ot), _pread(ob)]))
        self.assertEqual(counts.mod, 1)
        self.assertEqual(counts.nonmeth_unmod, 1)


# ---------------------------------------------------------------------------
# validate_bam_and_contig — fake BAM / FASTA handles
# ---------------------------------------------------------------------------

@dataclass
class FakeBam:
    """Mimics the subset of pysam.AlignmentFile used by validate_bam_and_contig."""
    _has_index: bool = True
    references: tuple = ("chrM",)

    def has_index(self) -> bool:
        return self._has_index


@dataclass
class FakeFasta:
    """Mimics the subset of pysam.FastaFile used by validate_bam_and_contig."""
    references: tuple = ("chrM",)


class TestValidateBamAndContig(unittest.TestCase):
    def _args(self, contig: str = "chrM") -> SimpleNamespace:
        return SimpleNamespace(
            contig=contig,
            bam=Path("/tmp/x.bam"),
            fasta=Path("/tmp/x.fa"),
        )

    def test_missing_index_exits(self):
        with self.assertRaises(SystemExit):
            smom.validate_bam_and_contig(FakeBam(_has_index=False), FakeFasta(), self._args())

    def test_contig_missing_from_bam_exits(self):
        bam = FakeBam(references=("chr1",))
        fasta = FakeFasta(references=("chrM",))
        with self.assertRaises(SystemExit):
            smom.validate_bam_and_contig(bam, fasta, self._args())

    def test_contig_missing_from_fasta_exits(self):
        bam = FakeBam(references=("chrM",))
        fasta = FakeFasta(references=("chr1",))
        with self.assertRaises(SystemExit):
            smom.validate_bam_and_contig(bam, fasta, self._args())

    def test_all_present_no_exit(self):
        # Should not raise.
        smom.validate_bam_and_contig(FakeBam(), FakeFasta(), self._args())


# ---------------------------------------------------------------------------
# Schema constancy — guards against accidental column drift
# ---------------------------------------------------------------------------

class TestSchemaConstancy(unittest.TestCase):
    def test_position_tsv_header_length(self):
        self.assertEqual(len(smom.POSITION_TSV_HEADER), 19)

    def test_position_tsv_last_column_is_beta(self):
        self.assertEqual(smom.POSITION_TSV_HEADER[-1], "beta_mod_fraction")

    def test_position_tsv_contains_required_columns(self):
        required = {
            "chrom",
            "pos_1based",
            "group",
            "ref_base",
            "mod_base",
            "unmod_base",
            "mod",
            "unmod",
            "other",
            "mod_fraction",
            "nonmeth_mod",
            "nonmeth_unmod",
            "nonmeth_mod_fraction",
            "beta_mod_fraction",
        }
        self.assertTrue(required.issubset(set(smom.POSITION_TSV_HEADER)))


if __name__ == "__main__":
    unittest.main(verbosity=2)
