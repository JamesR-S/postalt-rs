# Validation

## Reference comparison

The reference implementation is `bwa-postalt.js` from [BWA/bwakit](https://github.com/lh3/bwa/blob/master/bwakit/bwa-postalt.js). Verify that the script and k8 runtime are compatible before comparing outputs: older JavaScript byte-array operations can behave differently in newer k8 releases, especially for reverse-complemented supplementary alignments.

The required comparison is byte-for-byte identity of the complete SAM output, not just matching record counts or CIGAR strings.

```bash
python3 scripts/compare.py \
  --k8 /path/to/k8 \
  --js /path/to/compatible-bwa-postalt.js \
  --rust target/release/postalt-rs \
  --threads 4 \
  --alt reference.fa.alt \
  --sam pre_postalt.sam
```

`pre_postalt.sam` should contain SAM headers and complete, adjacent QNAME/mate-segment groups. Generate it directly from BWA before postalt. Do not create a test subset by cutting in the middle of a group.

## Parallel comparison

The parallel implementation must produce the same bytes as `-t 1` for each worker count, including output records emitted as supplementary alignments.

```bash
bash scripts/compare_parallel.sh \
  target/release/postalt-rs \
  reference.fa.alt \
  pre_postalt.sam \
  results/parallel
```

This checks `-t 2`, `-t 4` and `-t 8` against `-t 1`.

## Test coverage

The included Rust unit tests cover coordinate lifting, CIGAR handling, reverse complements, and serial/parallel output comparisons over several batch boundaries. Additional validation should include:

- ALT-associated and ordinary alignments in the same stream.
- Reverse-strand supplementary records with nonuniform quality strings.
- `XA:Z`, `lt:Z`, and `om:i` tag behaviour, including multiple hits.
- Insertions, deletions, clipping and HLA-region reads.
- Records spanning parallel batch boundaries.
- Missing or unusual optional SAM tags.
- Representative GRCh38+ALT WGS reads, including a complete large sample.

A reverse-complement test with identical quality characters cannot detect incorrect QUAL reversal; use varied quality characters for that check.

## Results available so far

Four synthetic datasets were compared against a k8-compatible JavaScript reference using the serial Rust implementation, with byte-identical SAM outputs:

| Dataset | Output SAM records |
| --- | ---: |
| Edge-case fixture | 12 |
| Smoke | 506 |
| Mixed | 9,361 |
| ALT-enriched | 12,410 |


## Postalt and mate tags

The implementation does not synchronize `MQ:i` tags across mates after changing MAPQ. For paired-end pipelines that retain mate mapping-quality tags, run `samtools fixmate` on name-sorted or name-collated postalt records. Check that downstream steps do not discard records unexpectedly.
