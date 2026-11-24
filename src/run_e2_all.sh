#!/bin/bash
# Engine 2 (Groebner/saturation) over all affine survivors, resumable: skips chunks with a .done marker.
cd "$(dirname "$0")/.."
for f in data/E2_chunks/in_*.txt; do
  b=${f%.txt}; o=${b/in_/out_}.jsonl
  [ -f "$o.done" ] && continue
  ALGEBRA_MAX_WORKERS=4 .venv/Scripts/python src/algebra_batch.py "$f" "$o" --engines 2 --workers 4 --timeout 120 && touch "$o.done"
  echo "$(date) $o" >> data/E2_chunks/progress.log
done
