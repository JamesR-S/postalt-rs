# Benchmarking

Benchmark post-processing on the **same pre-aligned SAM input** for both implementations. Sending output to `/dev/null` isolates most postalt work from downstream sort or compression costs.

```bash
bash scripts/benchmark.sh \
  /path/to/k8 \
  /path/to/bwa-postalt.js \
  reference.fa.alt \
  pre_postalt.sam \
  target/release/postalt-rs
```

This script defaults to serial Rust processing. To compare Rust worker counts:

```bash
bash scripts/compare_parallel.sh \
  target/release/postalt-rs \
  reference.fa.alt \
  pre_postalt.sam \
  results/parallel
```

The comparison script checks output equality before reporting timings.

## Serial reference measurement

One synthetic workload of 100,000 paired reads, measured over five runs, produced the following median postalt-only times:

| Implementation | Median wall time |
| --- | ---: |
| k8 + compatible `bwa-postalt.js` | 1.700 s |
| Serial Rust implementation | 0.850 s |

This is a **2.00× speedup on that workload**. It was measured with the earlier serial implementation, not the parallel release. It is not a GRCh38 throughput estimate.

