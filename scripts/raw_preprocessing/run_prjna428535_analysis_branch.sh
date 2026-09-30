#!/usr/bin/env bash
set -euo pipefail
ROOT="${STAGE1_ROOT:?Set STAGE1_ROOT}"
while [[ ! -s "$ROOT/logs/PRJNA428535_DADA2_COMPLETE" ]]; do
  if [[ -s "$ROOT/logs/PRJNA428535_DADA2_FAILED" ]]; then
    echo "DADA2 failed; not starting paired analysis" >&2
    exit 1
  fi
  sleep 30
done
exec python3 "$ROOT/workflow/analyze_prjna428535.py"
