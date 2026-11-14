#!/usr/bin/env bash
# Usage: compare_parallel.sh POSTALT_RS REF.fa.alt input.sam [OUTPUT_PREFIX]
# The script compares serial output with 2, 4, and 8 parallel workers.
set -euo pipefail
if [[ $# -lt 3 || $# -gt 4 ]]; then
  echo "Usage: $0 /path/postalt-rs reference.fa.alt input.sam [output-prefix]" >&2
  exit 2
fi
RUST="$1"
ALT="$2"
SAM="$3"
PREFIX="${4:-parallel-check}"
mkdir -p -- "$(dirname -- "$PREFIX")"
"$RUST" -t 1 "$ALT" "$SAM" > "${PREFIX}.serial.sam"
for n in 2 4 8; do
  echo "[check] $n workers: byte-for-byte equivalence..."
  "$RUST" -t "$n" "$ALT" "$SAM" > "${PREFIX}.t${n}.sam"
  if ! cmp -s "${PREFIX}.serial.sam" "${PREFIX}.t${n}.sam"; then
    echo "FAILED with $n workers; first differing lines:" >&2
    diff -u "${PREFIX}.serial.sam" "${PREFIX}.t${n}.sam" | head -80 || true
    exit 1
  fi
  echo "PASS: $n workers"
done

echo "[timing] postalt only, output discarded"
for n in 1 2 4 8; do
  for rep in 1 2 3; do
    /usr/bin/time -f "workers=$n run=$rep wall=%e rss_kb=%M" \
      "$RUST" -t "$n" "$ALT" "$SAM" > /dev/null
  done
done
