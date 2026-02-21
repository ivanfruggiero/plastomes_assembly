#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# 03_getorganelle.sh
# Run GetOrganelle (plastid assembly) on paired trimmed reads.
#
# Input naming convention (trimmed paired):
#   results/02_trimming/paired/SAMPLE_R1_paired.fq.gz
#   results/02_trimming/paired/SAMPLE_R2_paired.fq.gz
#
# Output:
#   results/03_assembly/getorganelle/SAMPLE/   (GetOrganelle output folder)
# Logs:
#   logs/03_getorganelle/SAMPLE.log
#
# Usage:
#   bash scripts/03_getorganelle.sh
#
# Optional:
#   --config  path/to/config.env     (default: config/config.env)
#   --threads N                      (overrides THREADS in config)
#   --conda                          (attempt to activate a conda env)
#   --conda-env NAME                 (default: GETORGANELLE_ENV_NAME in config, or "getorganelle")
#   --exclude "SAMPLE1,SAMPLE2"      (comma-separated sample IDs)
#
# Notes:
# - Parameters (rounds, k-mer set, mode) should be adapted to your organism/data.
# - This script is designed to be portable and avoids hardcoded system paths.
# ============================================================

CONFIG="config/config.env"
OVERRIDE_THREADS=""
USE_CONDA="false"
CONDA_ENV=""
EXCLUDE_CSV=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --config) CONFIG="${2:-}"; shift 2;;
    --threads) OVERRIDE_THREADS="${2:-}"; shift 2;;
    --conda) USE_CONDA="true"; shift 1;;
    --conda-env) CONDA_ENV="${2:-}"; shift 2;;
    --exclude) EXCLUDE_CSV="${2:-}"; shift 2;;
    -h|--help) sed -n '1,160p' "$0"; exit 0;;
    *) echo "[ERROR] Unknown argument: $1" >&2; exit 1;;
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
THREADS="${THREADS:-8}"

GETORGANELLE_BIN="${GETORGANELLE_BIN:-get_organelle_from_reads.py}"
GETORGANELLE_MODE="${GETORGANELLE_MODE:-embplant_pt}"
ROUNDS="${ROUNDS:-10}"
KSET="${KSET:-21,45,65,85,105}"
GETORGANELLE_ENV_NAME="${GETORGANELLE_ENV_NAME:-getorganelle}"

# Override threads from CLI if provided
if [[ -n "$OVERRIDE_THREADS" ]]; then
  THREADS="$OVERRIDE_THREADS"
fi

# Resolve conda env name
if [[ -n "$CONDA_ENV" ]]; then
  GETORGANELLE_ENV_NAME="$CONDA_ENV"
fi

RESULTS_DIR="$ROOT_DIR/$RESULTS_DIR"
LOG_DIR="$ROOT_DIR/$LOG_DIR"

READS_DIR="$RESULTS_DIR/02_trimming/paired"
OUT_ROOT="$RESULTS_DIR/03_assembly/getorganelle"
GO_LOG_DIR="$LOG_DIR/03_getorganelle"

mkdir -p "$OUT_ROOT" "$GO_LOG_DIR"

# -------------------------
# Optional conda activation
# -------------------------
if [[ "$USE_CONDA" == "true" ]]; then
  if command -v conda >/dev/null 2>&1; then
    # Enable conda in non-interactive shells (best-effort)
    CONDA_BASE="$(conda info --base)"
    # shellcheck disable=SC1090
    source "$CONDA_BASE/etc/profile.d/conda.sh"
    conda activate "$GETORGANELLE_ENV_NAME" || {
      echo "[ERROR] Failed to activate conda env: $GETORGANELLE_ENV_NAME" >&2
      exit 1
    }
  else
    echo "[ERROR] --conda was requested but 'conda' is not available in PATH" >&2
    exit 1
  fi
fi

# -------------------------
# Sanity checks
# -------------------------
[[ -d "$READS_DIR" ]] || { echo "[ERROR] Reads dir not found: $READS_DIR" >&2; exit 1; }

if ! command -v "$GETORGANELLE_BIN" >/dev/null 2>&1; then
  echo "[ERROR] GetOrganelle not found in PATH: $GETORGANELLE_BIN" >&2
  echo "Set GETORGANELLE_BIN in config/config.env (or install GetOrganelle in your environment)." >&2
  exit 1
fi

# -------------------------
# Exclude list (optional)
# -------------------------
declare -A EXCLUDE=()
if [[ -n "$EXCLUDE_CSV" ]]; then
  IFS=',' read -r -a ex_arr <<< "$EXCLUDE_CSV"
  for s in "${ex_arr[@]}"; do
    s="$(echo "$s" | xargs)"
    [[ -n "$s" ]] && EXCLUDE["$s"]=1
  done
fi

# -------------------------
# Discover samples
# -------------------------
shopt -s nullglob
R1S=( "$READS_DIR"/*_R1_paired.fq.gz )

if [[ ${#R1S[@]} -eq 0 ]]; then
  echo "[ERROR] No R1 paired reads found in $READS_DIR (*_R1_paired.fq.gz)" >&2
  exit 1
fi

echo "========================================="
echo "[INFO] plastomes_assembly - GetOrganelle batch started"
echo "[INFO] Reads dir : $READS_DIR"
echo "[INFO] Output    : $OUT_ROOT"
echo "[INFO] Logs      : $GO_LOG_DIR"
echo "[INFO] Threads   : $THREADS"
echo "[INFO] Mode      : $GETORGANELLE_MODE"
echo "[INFO] Rounds    : $ROUNDS"
echo "[INFO] Kset      : $KSET"
echo "[INFO] Date      : $(date)"
echo "========================================="

# -------------------------
# Main loop
# -------------------------
for r1 in "${R1S[@]}"; do
  sample="$(basename "$r1" _R1_paired.fq.gz)"
  r2="$READS_DIR/${sample}_R2_paired.fq.gz"

  if [[ ${EXCLUDE["$sample"]+yes} == "yes" ]]; then
    echo "[WARN] Excluding sample: $sample"
    continue
  fi

  if [[ ! -f "$r2" ]]; then
    echo "[WARN] Missing R2 for $sample -> skipping"
    continue
  fi

  outdir="$OUT_ROOT/$sample"
  logfile="$GO_LOG_DIR/${sample}.log"
  mkdir -p "$outdir"

  # Resume logic: skip if a final assembly is already present
  if [[ -s "$outdir/final_assembly_graph.fastg" ]] || compgen -G "$outdir/*.fasta" >/dev/null; then
    echo "[SKIP] $sample (output already present)"
    continue
  fi

  echo "[RUN]  $sample"
  echo "  R1: $r1"
  echo "  R2: $r2"
  echo "  OUT: $outdir"
  echo "  LOG: $logfile"

  "$GETORGANELLE_BIN" \
    -1 "$r1" -2 "$r2" \
    -o "$outdir" \
    -F "$GETORGANELLE_MODE" \
    -R "$ROUNDS" \
    -k "$KSET" \
    -t "$THREADS" \
    >"$logfile" 2>&1 || {
      echo "[FAIL] $sample (see $logfile)"
      continue
    }

  echo "[OK]   $sample"
done

echo "[DONE] GetOrganelle batch finished"
