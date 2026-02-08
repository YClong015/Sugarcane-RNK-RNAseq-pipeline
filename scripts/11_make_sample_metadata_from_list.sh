#!/usr/bin/env bash
set -euo pipefail
source scripts/utils.sh
load_cfg

[[ -f "$SAMPLES_RKN_ONLY_TXT" ]] || die "Missing: $SAMPLES_RKN_ONLY_TXT"

: "${CONTROL_CODE:?Missing CONTROL_CODE in config.env}"
: "${CONTROL_LABEL:?Missing CONTROL_LABEL in config.env}"
: "${CASE_CODE:?Missing CASE_CODE in config.env}"
: "${CASE_LABEL:?Missing CASE_LABEL in config.env}"

OUT_CSV="${DESEQ_DIR}/sample_metadata.csv"
echo "sample_id,genotype,time_point,treatment,replicate" > "$OUT_CSV"

awk -F'_' -v cc="$CONTROL_CODE" -v cl="$CONTROL_LABEL" \
  -v kc="$CASE_CODE" -v kl="$CASE_LABEL" 'BEGIN{OFS=","}
{
  rep=$4;
  for(i=5;i<=NF;i++) rep=rep"_"$i;

  tr=$3;
  if (tr == cc) tr = cl;
  if (tr == kc) tr = kl;

  print $0,$1,$2,tr,rep
}' "$SAMPLES_RKN_ONLY_TXT" >> "$OUT_CSV"

head "$OUT_CSV" || true
note "Metadata written: $OUT_CSV"
