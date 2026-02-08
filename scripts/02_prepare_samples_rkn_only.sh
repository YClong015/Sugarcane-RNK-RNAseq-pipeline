#!/usr/bin/env bash
set -euo pipefail
source scripts/utils.sh
load_cfg

# If missing, copy the repo template into $DESEQ_DIR and stop.
if [[ ! -f "$SAMPLES_RKN_ONLY_TXT" ]]; then
  tmpl="samples/samples_rkn_only.txt"
  [[ -f "$tmpl" ]] || die "Missing template: $tmpl"
  ensure_dir "$DESEQ_DIR"
  cp "$tmpl" "$SAMPLES_RKN_ONLY_TXT"
  note "Created: $SAMPLES_RKN_ONLY_TXT"
  note "Edit it, then rerun this script."
  exit 0
fi

# Sanitize: drop RLN lines, comments, empty lines; keep first column only.
tmp="${SAMPLES_RKN_ONLY_TXT}.tmp"
grep -v -E 'RLN|SES208_12w_RLN_3' "$SAMPLES_RKN_ONLY_TXT" \
  | grep -v -E '^[[:space:]]*#' \
  | awk 'NF > 0 { print $1 }' > "$tmp"

mv "$tmp" "$SAMPLES_RKN_ONLY_TXT"

if grep -q -E 'RLN|SES208_12w_RLN_3' "$SAMPLES_RKN_ONLY_TXT"; then
  die "RLN still present in $SAMPLES_RKN_ONLY_TXT"
fi

note "Sanitized: $SAMPLES_RKN_ONLY_TXT"
wc -l "$SAMPLES_RKN_ONLY_TXT" || true
