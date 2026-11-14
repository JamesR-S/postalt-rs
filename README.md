# postalt-rs

`postalt-rs` is a Rust implementation of [bwakit `bwa-postalt.js`](https://github.com/lh3/bwa/blob/master/bwakit/bwa-postalt.js) for post-processing BWA alignments against references containing ALT contigs.

It reads SAM from standard input (or a SAM file) and writes processed SAM to standard output. ALT-to-primary mappings are loaded from the reference `.alt` file. The implementation has no third-party Rust dependencies.

## Build

Rust 1.63 or later is required.

```bash
cargo test --offline
cargo build --release --offline
```

The executable is written to `target/release/postalt-rs`. A GitHub Actions workflow also builds a static Linux x86-64 executable.

## Usage

```text
postalt-rs [-t THREADS] [-r MIN_PA_RATIO] [-p HLA_PREFIX] reference.fa.alt [input.sam]
```

Without `input.sam`, alignments are read from standard input. Output is always SAM text on standard output.

| Option | Description |
| --- | --- |
| `-t N`, `--threads N` | Number of processing workers, 1–64 (default: 1) |
| `-r FLOAT` | Primary/ALT mapping-quality ratio threshold (default: 1) |
| `-p PREFIX` | Write HLA FASTQ files; processing runs in serial mode |
| `-h`, `--help` | Show command usage |
| `-v`, `--version` | Print the upstream-compatible `r985` marker |

For example:

```bash
postalt-rs -t 4 reference.fa.alt < input.sam > processed.sam
```

With bwa-mem3:

```bash
set -euo pipefail

bwa-mem3 mem -t 8 -M reference.fa R1.fastq.gz R2.fastq.gz \
  | postalt-rs -t 4 reference.fa.alt \
  | samtools sort -n -@ 3 -l 1 -O BAM -o aligned_ns.bam -
```

Input must contain complete, adjacent alignment groups for each read name and mate segment, as emitted by BWA. Arbitrarily reordered or coordinate-sorted SAM input is not supported. Each group is processed as a unit; parallel workers share the immutable ALT mapping and output batches are written in input order.

### Mate mapping quality

ALT post-processing can change an alignment's MAPQ without updating the `MQ:i` tag on its mate. This implementation preserves upstream postalt tag behaviour. If the aligner emits `MQ:i`, run `samtools fixmate` **after** postalt and after arranging the records in name order. `fixmate` updates mate-related tags; it does not recompute MAPQ itself.

### HLA extraction

The `-p` option writes per-HLA FASTQ files. Even if `-t` is greater than 1, HLA extraction uses the serial processing path.

## Tests

Run the Rust unit tests:

```bash
cargo test --offline
```

Compare the serial and parallel modes against a complete pre-postalt SAM:

```bash
bash scripts/compare_parallel.sh \
  ./target/release/postalt-rs \
  reference.fa.alt \
  input.sam \
  results/parallel
```

Compare with a k8-compatible copy of the original JavaScript implementation:

```bash
python3 scripts/compare.py \
  --k8 /path/to/k8 \
  --js /path/to/bwa-postalt.js \
  --rust ./target/release/postalt-rs \
  --threads 4 \
  --alt reference.fa.alt \
  --sam input.sam
```

The comparison is byte for byte. See [Validation](docs/VALIDATION.md) for datasets and checks, and [Benchmarking](docs/BENCHMARKING.md) for timing methods.

## Design

- SAM streaming with buffered input and output.
- Read-name and mate-segment grouping before processing.
- A read-only interval index for ALT mappings, shared across workers.
- Bounded processing batches and ordered output in parallel mode.
- Serial mode retained as the reference execution path.
- No native BAM input: convert BAM to SAM with `samtools view -h` if required.

## License and attribution

Licensed under GPL-3.0-or-later. See [COPYING](COPYING).

The underlying ALT post-processing algorithm originates in Heng Li's [BWA/bwakit](https://github.com/lh3/bwa) project. This is a separate Rust implementation and is not an official BWA release.
