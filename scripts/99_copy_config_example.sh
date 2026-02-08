#!/usr/bin/env bash
set -euo pipefail

src="config/config.env.example"
dst="config/config.env"

if [[ -f "$dst" ]]; then
  echo "[ERROR] $dst already exists." >&2
  exit 1
fi

if [[ ! -f "$src" ]]; then
  echo "[ERROR] Missing $src" >&2
  exit 1
fi

cp "$src" "$dst"
echo "[INFO] Created $dst"
