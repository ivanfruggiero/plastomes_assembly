#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# 01_qc.sh
# Run FastQC on paired-end reads (pre- or post-trimming).
#
# Naming conventions assumed:
#   - pre-trimming (raw):   SAMPLE_R1_001.fastq.gz / SAMPLE_R2_001.fastq.gz
#   - post-trimming:        SAMPLE_R1_paired.fq.gz / SAMPLE_R2_paired.fq.gz
#
# Usage:
#   bash scripts/01_qc.sh --stage pre
#   bash scripts/01_qc.sh --stage post
#
# Optional:
#   --config  path/to/config.env   (default: config/config.env)
#   --threads N                    (overrides THREADS in config)
# ============================================================

STAGE="pre"
CONFIG="config/config.env"
OVERRIDE_THREADS=""

# -------------------------
# Parse arguments
# -------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --stage)
      STAGE="${2:-}"
      shift 2
      ;;
    --config)
      CONFIG="${2:-}"
      shift 2
      ;;
    --threads)
      OVERRIDE_THREADS="${2:-}"
      shift 2
      ;;
    -h|--help)
      sed -n '1,80p' "$0"
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

RESULTS_DIR="${RESULTS_DIR:-results}"
LOG_DIR="${LOG_DIR:-logs}"
RAW_DIR="${RAW_DIR:-data/raw}"
FASTQC_BIN="${FASTQC_BIN:-fastqc}"
THREADS="${THREADS:-4}"

# Override threads from CLI if provided
if [[ -n "$OVERRIDE_THREADS" ]]; then
  THREADS="$OVERRIDE_THREADS"
fi

RESULTS_DIR="$ROOT_DIR/$RESULTS_DIR"
LOG_DIR="$ROOT_DIR/$LOG_DIR"
RAW_DIR="$ROOT_DIR/$RAW_DIR"

mkdir -p "$RESULTS_DIR" "$LOG_DIR"

# -------------------------
# Stage-specific settings
# -------------------------
case "$STAGE" in
  pre)
    IN_DIR="$RAW_DIR"
    OUT_DIR="$RESULTS_DIR/01_qc/pre_fastqc"
    LOGFILE="$LOG_DIR/01_qc_pre.log"
    R1_PATTERN="*_R1_001.fastq.gz"
    ;;
  post)
    IN_DIR="$RESULTS_DIR/02_trimming/paired"
    OUT_DIR="$RESULTS_DIR/01_qc/post_fastqc"
    LOGFILE="$LOG_DIR/01_qc_post.log"
    R1_PATTERN="*_R1_paired.fq.gz"
    ;;
  *)
    echo "[ERROR] --stage must be 'pre' or 'post' (got: $STAGE)" >&2
    exit 1
    ;;
esac

mkdir -p "$OUT_DIR"

# -------------------------
# Sanity checks
# -------------------------
if ! command -v "$FASTQC_BIN" >/dev/null 2>&1; then
  echo "[ERROR] FastQC not found in PATH: $FASTQC_BIN" >&2
  echo "Set FASTQC_BIN in config/config.env (e.g., FASTQC_BIN=/path/to/fastqc)" >&2
  exit 1
fi

if [[ ! -d "$IN_DIR" ]]; then
  echo "[ERROR] Input directory not found: $IN_DIR" >&2
  exit 1
fi

{
  echo "========================================="
  echo "[INFO] plastomes_assembly - QC started"
  echo "[INFO] Stage      : $STAGE"
  echo "[INFO] Input dir  : $IN_DIR"
  echo "[INFO] Output dir : $OUT_DIR"
  echo "[INFO] FastQC     : $FASTQC_BIN"
  echo "[INFO] Threads    : $THREADS"
  echo "[INFO] Date       : $(date)"
  echo "========================================="
} | tee "$LOGFILE"

# -------------------------
# Discover samples from R1 files
# -------------------------
shopt -s nullglob
R1_FILES=( "$IN_DIR"/$R1_PATTERN )

if [[ ${#R1_FILES[@]} -eq 0 ]]; then
  echo "[ERROR] No R1 files found in $IN_DIR matching: $R1_PATTERN" | tee -a "$LOGFILE"
  exit 1
fi

# -------------------------
# Main loop
# -------------------------
for r1 in "${R1_FILES[@]}"; do
  base="$(basename "$r1")"

  # Derive sample name and expected R2 based on stage naming convention
  if [[ "$STAGE" == "pre" ]]; then
    sample="${base%_R1_001.fastq.gz}"
    r2="$IN_DIR/${sample}_R2_001.fastq.gz"
  else
    sample="${base%_R1_paired.fq.gz}"
    r2="$IN_DIR/${sample}_R2_paired.fq.gz"
  fi

  if [[ ! -f "$r2" ]]; then
    echo "[WARN] Missing R2 for sample '$sample' -> skipping" | tee -a "$LOGFILE"
    continue
  fi

  echo "[INFO] Running FastQC for: $sample" | tee -a "$LOGFILE"

  "$FASTQC_BIN" \
    --threads "$THREADS" \
    --outdir "$OUT_DIR" \
    "$r1" "$r2" >>"$LOGFILE" 2>&1

done

echo "=========================================" | tee -a "$LOGFILE"
echo "[INFO] QC completed successfully" | tee -a "$LOGFILE"
echo "[INFO] End time: $(date)" | tee -a "$LOGFILE"
echo "=========================================" | tee -a "$LOGFILE"
