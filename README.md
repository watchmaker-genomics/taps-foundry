# taps-foundry

[rastair](https://bitbucket.org/bsblabludwig/rastair/src/master/) and [asTair](https://bitbucket.org/bsblabludwig/astair/src/master/) were created by Benjamin Schuster-Böckler's lab at the University of Oxford, Ludwig Institute for Cancer Research. 

Python scripts in this repository rely on Python 3+.

## `call-methylation.sh`

This script processes paired-end FASTQ files for TAPS methylation calling. It's intended as an example pipeline to help users get started with TAPS methylation analysis. Most users will likely pull out specific steps for their own custom pipelines, however it can be run as-is or with minimal modifications. See the **Usage** section below for details.

#### Usage

```
call-methylation.sh \
  --prefix SAMPLE \
  --cpus 16 \
  --memory 64 \
  --reference /path/to/reference.fasta \
  --adapter1 <seq> \
  --adapter2 <seq> \
  [--mapq 20] [--baseq 30]
```

Input files must be named:
```
<PREFIX>_R1_001.fastq.gz
<PREFIX>_R2_001.fastq.gz
```

#### Key Arguments

| Argument | Description |
|---------|-------------|
| `--prefix` | Sample name prefix |
| `--cpus` | Number of threads |
| `--memory` | RAM for GATK (GB) |
| `--reference` | Reference FASTA |
| `--adapter1`, `--adapter2` | Adapter sequences for trimming |

#### Variables you may wish to Adjust

| Variable | Script Value | Notes |
|----------|---------|-------|
| Trim length (`cutadapt -m`) | 20 | Minimum read length to retain through cutadapt |
| Read length | 151 | Adjust to fit your read length |
| asTair `start_clip` / `end_clip` | 7 / 7 | mbias mask for start/end of r1/r2 |
| rastair `nOT`/`nOB`| `0,0,20,0` | mbias mask for start/end of r1/r2 for OTs and OBs |

#### Example with Illumina Truseq Adapters

```
call-methylation.sh \
    --prefix=prefix \
    --cpus=40 \
    --memory=80 \
    --reference=/path/to/hg38.fa \
    --adapter1=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
    --adapter2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    --mapq=20 \
    --baseq=30
```

## `rastair_summarize.py`

This script summarizes methylation counts from a rastair `.mods` file.  
It filters out rows containing SNPs and outputs a one-line summary including:

- Total modified counts (`mod`)
- Total unmodified counts (`unmod`)
- Overall methylation rate (`mod / (mod + unmod)`)
- Number of covered CpG positions

#### Input Requirements
The input file must be a tab-delimited `.mods` file produced by `rastair call` and must contain the fields:
`#chr`, `start`, `end`, `name`, `beta_est`, `strand`, `unmod`, `mod`, `no_snp`, `snp`, `coverage`, `genotype`, `gt_p_score`, `gt_conf_score`.

#### Usage

```
python rastair_summarize.py \
    --input-file sample_rastair.mods \
    --output-dir results/
```

#### Output
Produces a `<input>.summary` file containing tab-separated fields:

```
total_mod    total_unmod    methylation_rate    covered_positions
```

Example:

```
total_mod    total_unmod    methylation_rate    covered_positions
15234        48791          0.23892             8462
```

## `save_as_methylkit.py`

This script converts a `rastair call` `.mods` file into a **methylKit-compatible** TSV format. It computes percent methylation (`freqC`) and percent unmethylated (`freqT`) per CpG position and outputs the standard methylKit columns. 

Note that with TAPS chemistry, `freqC` represents the count of mCtoT conversions which is the opposite of Bisulfite sequencing.

#### Input Requirements

The input file must be a tab-separated rastair mods table containing at least:
```
#chr, start, mod, unmod, coverage[, strand]
```

#### Output Format
The output file contains:
```
chrBase    chr    base    strand    coverage    freqC    freqT
```
- `freqC` = percent methylated (`mod / (mod + unmod) * 100`)
- `freqT` = percent unmethylated

### Usage

```
python save_as_methylkit.py \
    --input_file sample_rastair.mods \
    --output_file sample_methylkit.tsv
```

#### Example

```
python save_as_methylkit.py -i test.mods -o results/test.mods.methylkit
```

## `scan_mito_ot_ob_mod_unmod.py`

Scans a single contig of an aligned BAM (default `chrM`) for TAPS-style methylation calls and writes per-position and per-group summary TSVs.

Mechanically:

- At each reference **C**, OT (Original Top) strand reads are counted as **mod** (`T` = methylated C → T) vs **unmod** (`C` = unconverted).
- At each reference **G**, OB (Original Bottom) strand reads are counted as **mod** (`A`) vs **unmod** (`G`).
- At the same positions, reads of the *opposite* strand are tallied as a background-rate control — base calls of the same shape that can only come from sequencing error, deamination damage, or somatic SNV. Used to produce an error-corrected methylation estimate (`beta_mod_fraction`).

Trim masks (`--nOT` / `--nOB`) are applied in **read orientation**, matching rastair v1 semantics; reverse-strand reads have their 5'/3' trims swapped relative to the BAM SEQ. Overlapping paired-end mates are deduplicated to one fragment per position (higher-BQ call wins) *after* the quality / flag filters, so a filtered-out higher-BQ mate doesn't clobber the keep-able one. Mapping-quality and base-quality filtering is pushed into the htslib pileup engine.

#### Dependencies

- Python ≥ 3.9
- `pysam` (`pip install pysam`)
- An indexed BAM (`.bai`) and an indexed FASTA (`.fai`)

The script fails fast with a clean error if any of those are missing, if the BAM is unindexed, or if the requested contig isn't in the BAM / FASTA header.

#### Usage

```
python scan_mito_ot_ob_mod_unmod.py \
    --bam sample.bam \
    --fasta hg38.fa \
    --contig chrM \
    --out-tsv sample.positions.tsv \
    --summary-tsv sample.summary.tsv
```

#### Key Arguments

| Argument | Default | Description |
|---|---|---|
| `--bam` | (required) | Sorted, indexed BAM. |
| `--fasta` | (required) | Indexed FASTA. |
| `--contig` | `chrM` | Reference contig to scan. Must be present in both BAM and FASTA. |
| `--min-mapq` | `20` | Minimum mapping quality (applied at the pileup engine). |
| `--min-baseq` | `30` | Minimum base quality (applied at the pileup engine). |
| `--include-flags` | `3` | Require all bits (`0x1` paired + `0x2` proper-pair). Use `0` for single-end. |
| `--exclude-flags` | `3852` | Drop unmapped / mate-unmapped / secondary / vendor-QC-fail / duplicate / supplementary. |
| `--nOT` | `0,0,0,0` | OT trim mask `r1_5p,r1_3p,r2_5p,r2_3p` (rastair v1, read-orientation). |
| `--nOB` | `0,0,0,0` | OB trim mask, same shape. |
| `--stranded-read` | `R1` | Which mate carries the original DNA strand. Use `R2` for PBAT / swapped libraries. Single-end ignores this. |
| `--max-depth` | `200000` | pysam pileup cap. |
| `--out-tsv` | (required) | Per-position TSV (per-row schema below). |
| `--summary-tsv` | (required) | Aggregate summary TSV (`OT` / `OB` / `combined` rows). |

#### Per-position TSV columns

| Column | Description |
|---|---|
| `chrom`, `pos_1based` | Reference coordinates (1-based). |
| `group` | `OT` (reference is C) or `OB` (reference is G). |
| `ref_base`, `mod_base`, `unmod_base` | Reference / modified-call / unmodified-call base. For OT: `C` / `T` / `C`. For OB: `G` / `A` / `G`. |
| `mod`, `unmod`, `other` | Methylated-strand call counts after filtering and per-fragment dedup. `other` is anything that is neither the `mod_base` nor the `unmod_base`. |
| `informative` | `mod + unmod`. |
| `depth` | `informative + other`. |
| `mod_fraction` | `mod / informative` — raw methylation estimate, **uncorrected**. |
| `nonmeth_mod`, `nonmeth_unmod`, `nonmeth_other` | Same tally computed from **opposite-strand** reads at the same position (the control channel). |
| `nonmeth_informative`, `nonmeth_depth` | Derived totals from the control channel. |
| `nonmeth_mod_fraction` | The noise / background-error rate at this position. |
| `beta_mod_fraction` | **Error-corrected methylation estimate**: `max(0, (mod_fraction − nonmeth_mod_fraction) / (1 − nonmeth_mod_fraction))`. Clamped to ≥ 0. Falls back to `mod_fraction` when no non-meth reads exist or when `nonmeth_mod_fraction = 1` (homozygous-alt SNV). |

#### Summary TSV columns

One row each for `OT`, `OB`, and `combined`. Schema mirrors the per-position TSV with aggregates instead of per-position values:

| Column | Description |
|---|---|
| `group` | `OT`, `OB`, or `combined`. |
| `ref_base`, `mod_base`, `unmod_base` | Per-group labels (`combined` shows `C/G`, `T/A`, `C/G`). |
| `positions` | Total reference positions of this group on the scanned contig. |
| `covered_positions` | Positions with ≥ 1 passing meth-strand read. |
| `mod`, `unmod`, `other`, `informative`, `depth`, `mod_fraction` | Aggregate methylated-strand counts. |
| `nonmeth_mod`, `nonmeth_unmod`, `nonmeth_other`, `nonmeth_informative`, `nonmeth_depth`, `nonmeth_mod_fraction` | Aggregate control-channel counts. |
| `beta_mod_fraction` | Aggregate error-corrected estimate. |

#### Tests

The script ships with 74 unit tests (stdlib `unittest`, no third-party deps). Run from the repo root:

```
python3 -m unittest discover -s tests -v
```

## Variant Calling for TAPS Data
TAPS introduces predictable C→T conversions during library preparation, which can confound standard variant callers by mimicking SNPs. Although TAPS preserves DNA integrity far better than bisulfite or EM-seq while simultaneously minimizing base changes, these systematic conversions still require specialized handling.

For accurate SNP and methylation-aware variant detection, use the TAPS+ Variant Caller (TVC), specifically developed for TAPS data:

[TVC on GitHub](https://github.com/watchmaker-genomics/TVC)

# License and Use
This project is licensed under the [MIT License](./LICENSE).

- You are free to use, copy, modify, and distribute this software in both
  commercial and non-commercial projects, subject to the terms of the MIT License.
- Attribution is required: please retain the copyright notice and license
  in any redistributed versions or substantial portions of this software.
- All software is provided **“AS IS”**, without warranties or conditions of any kind.

## Important Notes
- **Research Use Only:** This code is provided for research and development purposes.
  It is **not validated for diagnostic or clinical use**, and Watchmaker Genomics
  makes no representations about regulatory compliance.
- **No Warranty:** Use of this code is at your own risk. Watchmaker Genomics shall
  not be liable for any claims or damages arising from its use.
- **Contributions:** If you contribute code, you agree that your contributions
  will be licensed under the same MIT License.

By making a contribution to this project, I certify that:

(a) The contribution was created by me and I have the right to submit it under the MIT License; or
(b) The contribution is based upon previous work that is covered by the MIT License and I have the right to submit it under the same license; or
(c) The contribution was provided to me by someone who certified (a) or (b).

I understand and agree that this project and the contribution are public and that a record of the contribution (including my name and email) is maintained indefinitely.
