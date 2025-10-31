#!/bin/bash
set -euo pipefail

# Parse named arguments (GNU getopt)
TEMP=$(getopt -o p:c:m:r:q:b: \
  --long prefix:,cpus:,memory:,reference:,mapq:,baseq:,adapter1:,adapter2: -- "$@")
eval set -- "$TEMP"

# Defaults
MAPQ=20
BASEQ=30
ADAPTER1=""
ADAPTER2=""

# Extract options
while true; do
  case "$1" in
  -p | --prefix)
    PREFIX=$2
    shift 2
    ;;
  -c | --cpus)
    CPUS=$2
    shift 2
    ;;
  -m | --memory)
    MEMORY=$2
    shift 2
    ;;
  -r | --reference)
    REFERENCE=$2
    shift 2
    ;;
  -q | --mapq)
    MAPQ=$2
    shift 2
    ;;
  -b | --baseq)
    BASEQ=$2
    shift 2
    ;;
  --adapter1)
    ADAPTER1=$2
    shift 2
    ;;
  --adapter2)
    ADAPTER2=$2
    shift 2
    ;;
  --)
    shift
    break
    ;;
  *)
    echo "Internal error!"
    exit 1
    ;;
  esac
done

# Required args check
if [ -z "${PREFIX:-}" ] || [ -z "${CPUS:-}" ] || [ -z "${MEMORY:-}" ] || [ -z "${REFERENCE:-}" ] ||
  [ -z "${ADAPTER1:-}" ] || [ -z "${ADAPTER2:-}" ]; then
  echo "Usage: $0 --prefix=<sample_prefix> --cpus=<num_cpus> --memory=<GB> --reference=<fasta> --adapter1 <seq> --adapter2 <seq> [--mapq=<int>] [--baseq=<int>]"
  exit 1
fi

echo "Initiating methylation calling with:"
echo "  prefix=${PREFIX}"
echo "  cpus=${CPUS}"
echo "  memory=${MEMORY}G"
echo "  reference=${REFERENCE}"
echo "  mapq=${MAPQ}"
echo "  baseq=${BASEQ}"
echo "  adapter1=${ADAPTER1}"
echo "  adapter2=${ADAPTER2}"

# Dirs
WORK_DIR="${PREFIX}/work"
RESULTS_DIR="${PREFIX}/results"
ASTAIR_DIR="${RESULTS_DIR}/asTair"
RASTAIR_DIR="${RESULTS_DIR}/rastair"
mkdir -p "${WORK_DIR}" "${ASTAIR_DIR}" "${RASTAIR_DIR}"

# Trim adapters
echo "==> Trimming adapters..."
time cutadapt \
  --cores "${CPUS}" \
  -a "${ADAPTER1}" \
  -A "${ADAPTER2}" \
  -o "${WORK_DIR}/${PREFIX}_R1_001_trimmed.fastq.gz" \
  -p "${WORK_DIR}/${PREFIX}_R2_001_trimmed.fastq.gz" \
  -m 20 \
  "${PREFIX}_R1_001.fastq.gz" \
  "${PREFIX}_R2_001.fastq.gz" \
  >"${WORK_DIR}/cutadapt.log"

# Align to reference
echo "==> Aligning to reference genome..."
STEP_AFFIX="trimmed"
OUTPUT_FILE="${WORK_DIR}/${PREFIX}_${STEP_AFFIX}"
time bwa mem \
  -t "${CPUS}" \
  "${REFERENCE}" \
  "${WORK_DIR}/${PREFIX}_R1_001_trimmed.fastq.gz" \
  "${WORK_DIR}/${PREFIX}_R2_001_trimmed.fastq.gz" |
  samtools view --threads "${CPUS}" -o "${OUTPUT_FILE}.bam" -

# Sort the bam
echo "==> Sorting the BAM file..."
STEP_AFFIX="sorted"
INPUT_FILE="${OUTPUT_FILE}.bam"
OUTPUT_FILE="${OUTPUT_FILE}_${STEP_AFFIX}"
time samtools sort \
  -o "${OUTPUT_FILE}.bam" \
  -O bam \
  -T "${OUTPUT_FILE}.bam.temp" \
  "${INPUT_FILE}"

# Mark duplicates
echo "==> Marking duplicates..."
STEP_AFFIX="marked"
INPUT_FILE="${OUTPUT_FILE}.bam"
OUTPUT_FILE="${OUTPUT_FILE}_${STEP_AFFIX}"
time gatk \
  --java-options "-Xmx${MEMORY}g" \
  MarkDuplicates \
  --INPUT "${INPUT_FILE}" \
  --OUTPUT "${OUTPUT_FILE}.bam" \
  --METRICS_FILE "${OUTPUT_FILE}.bam.metrics" \
  --TMP_DIR .

# Filter duplicates/secondary; keep proper pairs
echo "==> Filtering duplicates, secondary alignments and keeping only proper pairs..."
STEP_AFFIX="deduped_filtered"
INPUT_FILE="${OUTPUT_FILE}.bam"
OUTPUT_FILE="${OUTPUT_FILE}_${STEP_AFFIX}"
time samtools view \
  -b \
  -f 3 \
  -F 1024 \
  -F 2048 \
  -q "${MAPQ}" \
  -o "${OUTPUT_FILE}.bam" \
  "${INPUT_FILE}"

# Index
echo "==> Indexing BAM file..."
time samtools index "${OUTPUT_FILE}.bam"

MBIAS_BAM="${OUTPUT_FILE}"
METHYLATION_CALLS_BAM="${OUTPUT_FILE}"

# asTair mbias
echo "==> Running asTair mbias..."
INPUT_FILE="${MBIAS_BAM}.bam"
time astair mbias \
  --N_threads "${CPUS}" \
  --no_information 0 \
  --read_length 151 \
  -i "./${INPUT_FILE}" \
  -f "${REFERENCE}" \
  -d "${ASTAIR_DIR}"

# asTair call
echo "==> Calling methylation with asTair..."
STEP_AFFIX="mCtoT_CpG"
OUTPUT_FILE="${ASTAIR_DIR}/$(basename "${METHYLATION_CALLS_BAM}")_${STEP_AFFIX}"
time astair call \
  --N_threads "${CPUS}" \
  --no_information 0 \
  --skip_clip_overlap true \
  --start_clip 7 \
  --end_clip 7 \
  --context CpG \
  --max_depth 100000 \
  --minimum_base_quality "${BASEQ}" \
  --minimum_mapping_quality "${MAPQ}" \
  -i "./${METHYLATION_CALLS_BAM}.bam" \
  -f "${REFERENCE}" \
  -d "${ASTAIR_DIR}"

# rastair mbias
echo "==> Determining rastair mbias..."
STEP_AFFIX="rastair_mbias"
INPUT_FILE="${MBIAS_BAM}"
OUTPUT_FILE="${RASTAIR_DIR}/$(basename "${INPUT_FILE}")_${STEP_AFFIX}"
time rastair mbias \
  --threads "${CPUS}" \
  --fasta-file "${REFERENCE}" \
  "${INPUT_FILE}.bam" \
  >"${OUTPUT_FILE}.tsv"

# rastair call
echo "==> Calling methylation with rastair..."
STEP_AFFIX="rastair"
INPUT_FILE="${METHYLATION_CALLS_BAM}"
OUTPUT_FILE="${RASTAIR_DIR}/$(basename "${METHYLATION_CALLS_BAM}")_${STEP_AFFIX}"
time rastair call \
  --nOT 0,0,20,0 \
  --nOB 0,0,20,0 \
  --min-baseq "${BASEQ}" \
  --min-mapq "${MAPQ}" \
  --fasta-file "${REFERENCE}" \
  "${INPUT_FILE}.bam" \
  >"${OUTPUT_FILE}.mods"

# Summarize the rastair output
time rastair_summarize.py \
  --rastair_in "${OUTPUT_FILE}.mods" \
  --summary_out "${OUTPUT_FILE}.summary"
