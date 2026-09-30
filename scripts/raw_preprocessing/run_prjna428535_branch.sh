#!/usr/bin/env bash
set -euo pipefail

ROOT="${STAGE1_ROOT:?Set STAGE1_ROOT}"
while [[ ! -s "$ROOT/logs/PRJNA428535_DOWNLOAD_COMPLETE" ]]; do
  if [[ -s "$ROOT/logs/PRJNA428535_DOWNLOAD_FAILED" ]]; then
    echo "download branch failed; not starting DADA2" >&2
    exit 1
  fi
  sleep 30
done

if [[ -s "$ROOT/logs/PRJNA428535_DADA2_COMPLETE" ]]; then
  echo "DADA2 already complete"
  exit 0
fi

exec "${STAGE0_RSCRIPT:-Rscript}" \
  "$ROOT/workflow/run_prjna428535_dada2.R"
