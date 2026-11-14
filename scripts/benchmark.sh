#!/usr/bin/env bash
# Usage: ./scripts/benchmark.sh /path/to/k8 /path/to/bwa-postalt.js /path/to/alt.sam /path/to/input.sam [path/to/postalt-rs]
set -euo pipefail
if [[ $# -lt 4 || $# -gt 5 ]]; then
  echo "Usage: $0 k8 bwa-postalt.js alt.sam input.sam [postalt-rs]" >&2
  exit 2
fi
k8=$1 js=$2 alt=$3 sam=$4 rs=${5:-target/release/postalt-rs}
for i in 1 2 3; do
  echo "--- pass $i: k8 ---" >&2
  /usr/bin/time -f 'wall %e s | user %U s | peak_rss %M KiB' "$k8" "$js" "$alt" "$sam" >/dev/null
  echo "--- pass $i: rust ---" >&2
  /usr/bin/time -f 'wall %e s | user %U s | peak_rss %M KiB' "$rs" "$alt" "$sam" >/dev/null
done
