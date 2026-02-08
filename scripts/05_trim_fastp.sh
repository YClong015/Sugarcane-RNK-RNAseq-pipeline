#!/usr/bin/env bash
set -euo pipefail
source scripts/utils.sh
load_cfg

need_cmd fastp

[[ -f "$SAMPLES_RKN_ONLY_TXT" ]] || die "Missing: $SAMPLES_RKN_ONLY_TXT"
ensure_dir "$CLEAN_DIR"

# Parameters (defined in config.env)
: "${FASTP_TRIM_FRONT1:?Missing FASTP_TRIM_FRONT1 in config.env}"
: "${FASTP_TRIM_FRONT2:?Missing FASTP_TRIM_FRONT2 in config.env}"
: "${FASTP_CUT_WINDOW_SIZE:?Missing FASTP_CUT_WINDOW_SIZE in config.env}"
: "${FASTP_CUT_MEAN_QUALITY:?Missing FASTP_CUT_MEAN_QUALITY in config.env}"
: "${FASTP_LENGTH_REQUIRED:?Missing FASTP_LENGTH_REQUIRED in config.env}"

while read -r s; do
  note "fastp: $s"

  fastp \
    --in1 "${MERGED_DIR}/${s}_R1.fastq.gz" \
    --in2 "${MERGED_DIR}/${s}_R2.fastq.gz" \
    --out1 "${CLEAN_DIR}/${s}_R1.fastq.gz" \
    --out2 "${CLEAN_DIR}/${s}_R2.fastq.gz" \
    --detect_adapter_for_pe \
    --trim_front1 "$FASTP_TRIM_FRONT1" \
    --trim_front2 "$FASTP_TRIM_FRONT2" \
    --cut_tail \
    --cut_window_size "$FASTP_CUT_WINDOW_SIZE" \
    --cut_mean_quality "$FASTP_CUT_MEAN_QUALITY" \
    --length_required "$FASTP_LENGTH_REQUIRED" \
    --thread "$THREADS" \
    --html "${CLEAN_DIR}/${s}.fastp.html" \
    --json "${CLEAN_DIR}/${s}.fastp.json"
done < "$SAMPLES_RKN_ONLY_TXT"

note "fastp finished: $CLEAN_DIR"

