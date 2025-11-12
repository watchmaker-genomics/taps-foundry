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
