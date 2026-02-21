#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# 02_trimming.sh
# Paired-end trimming with Trimmomatic.
#
# IMPORTANT:
#   Trimming settings MUST be adapted based on FastQC results
#   (and, if available, MultiQC summaries).
#   This script provides a reasonable, generic starting point.
#
# Input naming convention (raw):
#   SAMPLE_R1_001.fastq.gz / SAMPLE_R2_001.fastq.gz
#
# Output naming convention (trimmed):
#   results/02_trimming/paired/   SAMPLE_R1_paired.fq.gz / SAMPLE_R2_paired.fq.gz
#   results/02_trimming/unpaired/ SAMPLE_R1_unpaired.fq.gz / SAMPLE_R2_unpaired.fq.gz
#
# Usage:
#   bash scripts/02_trimming.sh
#
# Optional:
#   --config  path/to/config.env   (default: config/config.env)
#   --threads N                    (overrides THREADS in config)
#   --exclude "SAMPLE1,SAMPLE2"    (comma-separated sample IDs)
#
# Config keys required:
#   RAW_DIR
#   RESULTS_DIR
#   LOG_DIR
#   TRIMMOMATIC_JAR
#   ADAPTERS_FA
# ============================================================

CONFIG="config/config.env"
OVERRIDE_THREADS=""
EXCLUDE_CSV=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --config)
      CONFIG="${2:-}"
      shift 2
      ;;
    --threads)
      OVERRIDE_THREADS="${2:-}"
      shift 2
      ;;
    --exclude)
      EXCLUDE_CSV="${2:-}"
      shift 2
      ;;
    -h|--help)
      sed -n '1,120p' "$0"
      exit 0
      ;;
    *)
      echo "[ERROR] Unknown argument: $1" >&2
      exit 1
      ;;
  esac
done

# -------------------------
# Load config
# -------------------------
if [[ ! -f "$CONFIG" ]]; then
  echo "[ERROR] Config file not found: $CONFIG" >&2
  echo "Create it by copying config/config.example.env -> config/config.env" >&2
  exit 1
fi

# shellcheck disable=SC1090
source "$CONFIG"

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

RAW_DIR="${RAW_DIR:-data/raw}"
RESULTS_DIR="${RESULTS_DIR:-results}"
LOG_DIR="${LOG_DIR:-logs}"
THREADS="${THREADS:-8}"

TRIMMOMATIC_JAR="${TRIMMOMATIC_JAR:-}"
ADAPTERS_FA="${ADAPTERS_FA:-}"

# Override threads from CLI if provided
if [[ -n "$OVERRIDE_THREADS" ]]; then
  THREADS="$OVERRIDE_THREADS"
fi

RAW_DIR="$ROOT_DIR/$RAW_DIR"
RESULTS_DIR="$ROOT_DIR/$RESULTS_DIR"
LOG_DIR="$ROOT_DIR/$LOG_DIR"

OUT_BASE="$RESULTS_DIR/02_trimming"
PAIRED_DIR="$OUT_BASE/paired"
UNPAIRED_DIR="$OUT_BASE/unpaired"
TRIM_LOG_DIR="$LOG_DIR/02_trimming"

mkdir -p "$PAIRED_DIR" "$UNPAIRED_DIR" "$TRIM_LOG_DIR"

# -------------------------
# Sanity checks
# -------------------------
[[ -d "$RAW_DIR" ]] || { echo "[ERROR] RAW_DIR not found: $RAW_DIR" >&2; exit 1; }

[[ -n "$TRIMMOMATIC_JAR" ]] || {
  echo "[ERROR] TRIMMOMATIC_JAR is not set in config/config.env" >&2
  echo "Example: TRIMMOMATIC_JAR=/path/to/Trimmomatic-0.39/trimmomatic-0.39.jar" >&2
  exit 1
}
[[ -f "$TRIMMOMATIC_JAR" ]] || { echo "[ERROR] Trimmomatic jar not found: $TRIMMOMATIC_JAR" >&2; exit 1; }

[[ -n "$ADAPTERS_FA" ]] || {
  echo "[ERROR] ADAPTERS_FA is not set in config/config.env" >&2
  echo "Example: ADAPTERS_FA=/path/to/adapters/NexteraPE-PE.fa" >&2
  exit 1
}
[[ -f "$ADAPTERS_FA" ]] || { echo "[ERROR] Adapter file not found: $ADAPTERS_FA" >&2; exit 1; }

command -v java >/dev/null 2>&1 || { echo "[ERROR] java not found in PATH" >&2; exit 1; }

# -------------------------
# Exclude list (optional)
# -------------------------
declare -A EXCLUDE=()
if [[ -n "$EXCLUDE_CSV" ]]; then
  IFS=',' read -r -a ex_arr <<< "$EXCLUDE_CSV"
  for s in "${ex_arr[@]}"; do
    s="$(echo "$s" | xargs)"  # trim spaces
    [[ -n "$s" ]] && EXCLUDE["$s"]=1
  done
fi

# -------------------------
# Discover samples
# -------------------------
shopt -s nullglob
R1_FILES=( "$RAW_DIR"/*_R1_001.fastq.gz )

if [[ ${#R1_FILES[@]} -eq 0 ]]; then
  echo "[ERROR] No R1 files found in $RAW_DIR matching *_R1_001.fastq.gz" >&2
  exit 1
fi

# -------------------------
# Generic trimming parameters
# (ADAPT to FastQC)
# -------------------------
# Notes:
# - ILLUMINACLIP removes adapter contamination.
# - SLIDINGWINDOW removes low-quality tails.
# - MINLEN removes too-short reads.
# - Additional steps like HEADCROP / LEADING / TRAILING may be needed depending on QC.
ILLUMINACLIP_OPTS="ILLUMINACLIP:${ADAPTERS_FA}:2:30:10"
SLIDINGWINDOW_OPTS="SLIDINGWINDOW:4:20"
MINLEN_OPTS="MINLEN:50"
# Optional (disabled by default):
# LEADING_OPTS="LEADING:3"
# TRAILING_OPTS="TRAILING:3"
# HEADCROP_OPTS="HEADCROP:0"

{
  echo "========================================="
  echo "[INFO] plastomes_assembly - Trimming started"
  echo "[INFO] Raw dir        : $RAW_DIR"
  echo "[INFO] Paired out     : $PAIRED_DIR"
  echo "[INFO] Unpaired out   : $UNPAIRED_DIR"
  echo "[INFO] Threads        : $THREADS"
  echo "[INFO] Trimmomatic jar: $TRIMMOMATIC_JAR"
  echo "[INFO] Adapters       : $ADAPTERS_FA"
  echo "[INFO] Date           : $(date)"
  echo "========================================="
} | tee "$TRIM_LOG_DIR/02_trimming.run.log"

# -------------------------
# Main loop
# -------------------------
for r1 in "${R1_FILES[@]}"; do
  sample="$(basename "$r1" _R1_001.fastq.gz)"
  r2="$RAW_DIR/${sample}_R2_001.fastq.gz"

  if [[ ${EXCLUDE["$sample"]+yes} == "yes" ]]; then
    echo "[WARN] Excluding sample: $sample" | tee -a "$TRIM_LOG_DIR/02_trimming.run.log"
    continue
  fi

  if [[ ! -f "$r2" ]]; then
    echo "[WARN] Missing R2 for $sample -> skipping" | tee -a "$TRIM_LOG_DIR/02_trimming.run.log"
    continue
  fi

  echo "[INFO] Processing: $sample" | tee -a "$TRIM_LOG_DIR/02_trimming.run.log"

  sample_log="$TRIM_LOG_DIR/${sample}.trimmomatic.log"

  # Trimmomatic PE mode
  java -jar "$TRIMMOMATIC_JAR" PE \
    -threads "$THREADS" \
    "$r1" "$r2" \
    "$PAIRED_DIR/${sample}_R1_paired.fq.gz" \
    "$UNPAIRED_DIR/${sample}_R1_unpaired.fq.gz" \
    "$PAIRED_DIR/${sample}_R2_paired.fq.gz" \
    "$UNPAIRED_DIR/${sample}_R2_unpaired.fq.gz" \
    "$ILLUMINACLIP_OPTS" \
    "$SLIDINGWINDOW_OPTS" \
    "$MINLEN_OPTS" \
    >"$sample_log" 2>&1 || {
      echo "[ERROR] Trimming failed for $sample (see $sample_log)" | tee -a "$TRIM_LOG_DIR/02_trimming.run.log"
      continue
    }

  echo "[INFO] Finished: $sample" | tee -a "$TRIM_LOG_DIR/02_trimming.run.log"
done

echo "=========================================" | tee -a "$TRIM_LOG_DIR/02_trimming.run.log"
echo "[INFO] Trimming completed" | tee -a "$TRIM_LOG_DIR/02_trimming.run.log"
echo "[INFO] End time: $(date)" | tee -a "$TRIM_LOG_DIR/02_trimming.run.log"
echo "=========================================" | tee -a "$TRIM_LOG_DIR/02_trimming.run.log"
