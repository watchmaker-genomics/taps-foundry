#!/usr/bin/env python3
"""
This tool takes output from rastair and converts it to a VCF file from the perspective of variant calling.

Usage:
    python rastair_to_vcf.py -i input.rastair -s sample_name -o output.vcf -f reference.fa.fai [-p score_filter]
"""

import os
import sys
import datetime
import argparse
from csv import DictReader
from typing import Tuple
import enum

import pysam

# Define some constants that are referenced multiple times
RASTAIR_CONTIG = "#chr"
RASTAIR_POS = "start"
RASTAIR_BETA_EST = "beta_est"
RASTAIR_STRAND = "strand"
RASTAIR_UNMOD = "unmod"
RASTAIR_MOD = "mod"
RASTAIR_NO_SNP = "no_snp"
RASTAIR_SNP = "snp"
RASTAIR_COVERAGE = "coverage"
RASTAIR_GENOTYPE = "genotype"
RASTAIR_GT_P_SCORE = "gt_p_score"

EXPECTED_COLUMNS = {
    RASTAIR_CONTIG,
    RASTAIR_POS,
    RASTAIR_BETA_EST,
    RASTAIR_STRAND,
    RASTAIR_UNMOD,
    RASTAIR_MOD,
    RASTAIR_NO_SNP,
    RASTAIR_SNP,
    RASTAIR_COVERAGE,
    RASTAIR_GENOTYPE,
    RASTAIR_GT_P_SCORE,
}

VCF_MS = "MS"
VCF_GT = "GT"
VCF_BT = "BT"
VCF_DP = "DP"
VCF_MF = "MF"
VCF_MO = "MO"
VCF_NMO = "NMO"
VCF_SAF = "SAF"
VCF_SAO = "SAO"
VCF_SDP = "SDP"


class Filter(enum.Enum):
    """
    Filters for the VCF file
    """

    NO_VARIANT = "MethylationSiteWithoutVariant"
    LOW_SCORE = "SNPBelowScoreThreshold"
    PASS = "PASS"


class RastairRecord:
    """
    This class represents a single record from the Rastair output file.
    It knows how to represent itself as a VCF record.

    :param contig: The contig name
    :param pos: The position of the record
    :param beta_est: The beta estimate for the methylation rate
    :param strand: The strand of the record
    :param unmod: The number of unmodified reads
    :param mod: The number of modified reads
    :param no_snp: The number of reads without a SNP
    :param snp: The number of reads with a SNP
    :param coverage: The total coverage at the position
    :param rastair_genotype: The genotype from Rastair
    :param gt_p_score: The genotype p-score from Rastair
    """

    def __init__(
        self,
        contig: str,
        pos: int,
        beta_est: float,
        strand: str,
        unmod: int,
        mod: int,
        no_snp: int,
        snp: int,
        coverage: int,
        rastair_genotype: str,
        gt_p_score: float,
        score_threshold: int,
    ):
        self.contig = contig
        self.pos = pos
        self.beta_est = beta_est
        self.strand = strand
        self.unmod = unmod
        self.mod = mod
        self.no_snp = no_snp
        self.snp = snp
        self.coverage = coverage
        self.rastair_genotype = rastair_genotype
        self.gt_p_score = gt_p_score
        self._ref = None
        self._alt = None
        self._vcf_genotype = None
        self._score_threshold = score_threshold

    @property
    def ref(self) -> str:
        """
        The reference allele for the record
        """
        if self._ref is not None:
            return self._ref

        self._ref = "C" if self.strand == "+" else "G"
        return self._ref

    @property
    def alt(self) -> str:
        """
        The alternate allele for the record. If its not a snp it's represented as a methylation event
        """
        if self._alt is not None:
            return self._alt
        for allele in self.rastair_genotype.split("/"):
            if allele != self.ref:
                self._alt = allele
                return allele
        if self.strand == "+":
            self._alt = "T"
            return self._alt
        self._alt = "A"
        return self._alt

    @property
    def vcf_genotype(self) -> Tuple:
        """
        The genotype as a tuple of integers for the VCF record
        """
        if self._vcf_genotype is not None:
            return self._vcf_genotype
        alleles = self.rastair_genotype.split("/")
        self._vcf_genotype = tuple(0 if allele == self.ref else 1 for allele in alleles)
        return self._vcf_genotype

    def get_vcf_score(self) -> float:
        """
        The quality score for the VCF record from the perspective of mutation
        """
        if not self.is_snp():
            return 0.0
        return self.gt_p_score

    def is_snp(self) -> bool:
        """
        Is this a legitimate SNP
        """
        return 1 in self.vcf_genotype

    def retrieve_filter(self) -> Filter:
        """
        Retrieve the filter for the record

        :return: The filter pass is included
        """
        if not self.is_snp():
            return Filter.NO_VARIANT
        if self.gt_p_score < self._score_threshold:
            return Filter.LOW_SCORE
        return Filter.PASS

    def construct_vcf_record(
        self, sample_name: str, vcf_out: pysam.VariantFile
    ) -> pysam.VariantRecord:
        """
        Construct a VCF record from the Rastair record

        :param sample_name: The name of the sample
        :param vcf_out: The VCF file object
        :return: The VCF record
        """
        filter = self.retrieve_filter().value
        record = vcf_out.new_record(
            contig=self.contig,
            start=self.pos,
            stop=self.pos + 1,
            alleles=(self.ref, self.alt),
            id=".",
            qual=self.get_vcf_score(),
            filter=filter,
        )

        # Add INFO fields
        record.info[VCF_MS] = self.strand

        # Add FORMAT fields for each sample
        record.samples[sample_name][VCF_GT] = self.vcf_genotype
        record.samples[sample_name][VCF_BT] = self.beta_est
        record.samples[sample_name][VCF_DP] = self.coverage
        record.samples[sample_name][VCF_MF] = (
            self.mod / (self.mod + self.unmod) if self.mod + self.unmod > 0.0 else 0.0
        )
        record.samples[sample_name][VCF_MO] = self.mod
        record.samples[sample_name][VCF_NMO] = self.unmod
        record.samples[sample_name][VCF_SAF] = (
            self.snp / (self.snp + self.no_snp) if self.snp + self.no_snp > 0.0 else 0.0
        )
        record.samples[sample_name][VCF_SAO] = self.snp
        record.samples[sample_name][VCF_SDP] = self.snp + self.no_snp
        return record


class RastairVcfWriter:
    """
    This class is responsible for converting a Rastair file to a VCF file
    """

    def __init__(
        self, output_file: str, sample_name: str, fasta_index: str, score_filter: int
    ):
        """
        :param output_file: The output VCF file
        :param sample_name: The name of the sample
        :param fasta_index: The fasta index file
        :param score_filter: All scores below this value will be filtered out
        """

        self._score_filter = score_filter
        self._sample_name = sample_name
        header = self._configure_header(fasta_index, sample_name)
        self._vcf_out = pysam.VariantFile(output_file, "w", header=header)
        self._is_closed = False

    def __enter__(self):
        """
        Context manager entry point
        """
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        """
        Context manager exit point
        """
        self.close()

    def close(self):
        """
        Close the VCF file
        """
        if not self._is_closed:
            self._vcf_out.close()
            self._is_closed = True

    def _configure_header(
        self, fasta_index: str, sample_name: str
    ) -> pysam.VariantHeader:
        """
        Configure the VCF header

        :param fasta_index: The fasta index file
        :param sample_name: The name of the sample
        """
        header = pysam.VariantHeader()
        header.add_meta("fileformat", "VCFv4.2")
        header.add_meta("fileDate", datetime.datetime.now().strftime("%Y%m%d"))
        header.add_meta("source", "rastair_to_vcf")
        header.add_meta("reference", os.path.basename(fasta_index).split(".")[0])

        header.add_sample(sample_name)

        with open(fasta_index) as fai:
            for line in fai:
                tokens = line.split("\t")
                if len(tokens) < 2:
                    continue
                contig = tokens[0]
                length = tokens[1]
                header.contigs.add(contig, length=int(length))

        header.add_meta(
            "FILTER",
            items=[
                ("ID", Filter.NO_VARIANT.value),
                (
                    "Description",
                    "This is a potential methylation site without a snp being called",
                ),
            ],
        )

        header.add_meta(
            "FILTER",
            items=[
                ("ID", Filter.LOW_SCORE.value),
                (
                    "Description",
                    "This snp is not homozygous ref but has a score below the threshold",
                ),
            ],
        )

        header.add_meta(
            "INFO",
            items=[
                ("ID", VCF_MS),
                ("Number", "1"),
                ("Type", "String"),
                ("Description", "Expected Methylation Strand"),
            ],
        )

        header.add_meta(
            "FORMAT",
            items=[
                ("ID", "GT"),
                ("Number", "1"),
                ("Type", "String"),
                ("Description", "Genotype"),
            ],
        )

        header.add_meta(
            "FORMAT",
            items=[
                ("ID", VCF_BT),
                ("Number", "1"),
                ("Type", "Float"),
                ("Description", "Beta Estimate for methylation rate"),
            ],
        )
        header.add_meta(
            "FORMAT",
            items=[
                ("ID", VCF_DP),
                ("Number", "1"),
                ("Type", "Integer"),
                ("Description", "Total Depth"),
            ],
        )
        header.add_meta(
            "FORMAT",
            items=[
                ("ID", VCF_MF),
                ("Number", "1"),
                ("Type", "Float"),
                ("Description", "Methylation Frequency"),
            ],
        )
        header.add_meta(
            "FORMAT",
            items=[
                ("ID", VCF_MO),
                ("Number", "1"),
                ("Type", "Integer"),
                ("Description", "Methylation Observations"),
            ],
        )
        header.add_meta(
            "FORMAT",
            items=[
                ("ID", VCF_NMO),
                ("Number", "1"),
                ("Type", "Integer"),
                ("Description", "Non Methylation Observations"),
            ],
        )
        header.add_meta(
            "FORMAT",
            items=[
                ("ID", VCF_SAF),
                ("Number", "1"),
                ("Type", "Float"),
                ("Description", "SNP Allele Frequency on non methylated strand"),
            ],
        )
        header.add_meta(
            "FORMAT",
            items=[
                ("ID", VCF_SAO),
                ("Number", "1"),
                ("Type", "Integer"),
                (
                    "Description",
                    "SNP Alternate Allele Observation on non methylated strand",
                ),
            ],
        )
        header.add_meta(
            "FORMAT",
            items=[
                ("ID", VCF_SDP),
                ("Number", "1"),
                ("Type", "Integer"),
                ("Description", "SNP depth on non methylated strand"),
            ],
        )
        return header

    def _iterate_records(self, rastair_file: str):
        """
        Iterate over the records in the Rastair file

        :param rastair_file: The Rastair file
        """
        with open(rastair_file) as f:
            reader = DictReader(f, delimiter="\t")
            for record in reader:
                if not EXPECTED_COLUMNS.issubset(record.keys()):
                    missed_columns = EXPECTED_COLUMNS - record.keys()
                    raise ValueError(
                        f"The header from provided Rastair file does not contain all expected columns. "
                        f"Here's what was missed: {missed_columns}"
                    )

                yield RastairRecord(
                    contig=record[RASTAIR_CONTIG],
                    pos=int(record[RASTAIR_POS]),
                    beta_est=float(record[RASTAIR_BETA_EST]),
                    strand=record[RASTAIR_STRAND],
                    unmod=int(record[RASTAIR_UNMOD]),
                    mod=int(record[RASTAIR_MOD]),
                    no_snp=int(record[RASTAIR_NO_SNP]),
                    snp=int(record[RASTAIR_SNP]),
                    coverage=int(record[RASTAIR_COVERAGE]),
                    rastair_genotype=record[RASTAIR_GENOTYPE],
                    gt_p_score=float(record[RASTAIR_GT_P_SCORE]),
                    score_threshold=self._score_filter,
                )

    def convert_file(self, rastair_file: str):
        """
        Convert the Rastair file to a VCF file

        :param rastair_file: The Rastair file
        """
        for record in self._iterate_records(rastair_file):
            vcf_record = record.construct_vcf_record(self._sample_name, self._vcf_out)
            self._vcf_out.write(vcf_record)


def parse_args():
    """
    Parse command line arguments
    """
    parser = argparse.ArgumentParser(description="Convert Rastair output to VCF format")
    parser.add_argument(
        "-i",
        "--input_rastair",
        required=True,
        help="The input rastair file",
    )
    parser.add_argument(
        "-s",
        "--sample_name",
        required=True,
        help="The sample name",
    )
    parser.add_argument(
        "-o",
        "--output_vcf",
        required=True,
        help="The output VCF file",
    )
    parser.add_argument(
        "-f",
        "--fasta_index",
        required=True,
        help="The fasta index file (.fai)",
    )
    parser.add_argument(
        "-p",
        "--score_filter",
        required=False,
        help="All scores below this value will be filtered out (default: 6)",
        type=int,
        default=6,
    )
    return parser.parse_args()


def main():
    """
    Main entry point for the standalone script
    """
    args = parse_args()
    with RastairVcfWriter(
        output_file=args.output_vcf,
        sample_name=args.sample_name,
        fasta_index=args.fasta_index,
        score_filter=args.score_filter,
    ) as writer:
        writer.convert_file(args.input_rastair)


if __name__ == "__main__":
    main()
