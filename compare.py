#!/usr/bin/env python3
"""Compare complete SAM outputs from original k8 bwa-postalt.js and postalt-rs.

Example:
 python3 scripts/compare.py --k8 /path/to/k8 --js /path/to/bwa-postalt.js \
   --rust target/release/postalt-rs --alt ref.fa.alt --sam input.sam
"""
import argparse
import os
import subprocess
import sys
import tempfile
from itertools import zip_longest
import time
from pathlib import Path


def execute(label, argv, stdout_path):
    t0 = time.perf_counter()
    with open(stdout_path, 'wb') as dest:
        result = subprocess.run(argv, stdout=dest, stderr=subprocess.PIPE)
    sec = time.perf_counter() - t0
    print(f'{label}: exit={result.returncode}, wall={sec:.2f}s, bytes={Path(stdout_path).stat().st_size:,}')
    if result.returncode:
        print(result.stderr.decode(errors='replace')[-4000:], file=sys.stderr)
        raise SystemExit(f'{label} failed')


def first_difference(a, b):
    with open(a, 'rb') as fa, open(b, 'rb') as fb:
        for i, (x, y) in enumerate(zip_longest(fa, fb, fillvalue=b''), 1):
            if x != y:
                return i, x.decode(errors='replace'), y.decode(errors='replace')
    return None


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--k8', required=True, help='Path to k8 binary')
    p.add_argument('--js', required=True, help='Path to original bwa-postalt.js')
    p.add_argument('--rust', default='target/release/postalt-rs')
    p.add_argument('--threads', type=int, default=1, help='Rust postalt worker threads (default: 1)')
    p.add_argument('--alt', required=True)
    p.add_argument('--sam', required=True)
    p.add_argument('--ratio', type=float, default=None)
    p.add_argument('--hla-prefix', default=None, help='Enable -p and compare HLA output files')
    args = p.parse_args()
    if not 1 <= args.threads <= 64:
        p.error('--threads must be between 1 and 64')
    js_opts = []
    rs_opts = []
    if args.ratio is not None:
        js_opts += ['-r', str(args.ratio)]
        rs_opts += ['-r', str(args.ratio)]
    with tempfile.TemporaryDirectory(prefix='postalt_diff_') as tmp:
        orig = os.path.join(tmp, 'original.sam')
        rust = os.path.join(tmp, 'rust.sam')
        if args.hla_prefix:
            js_opts += ['-p', os.path.join(tmp, 'js')]
            rs_opts += ['-p', os.path.join(tmp, 'rs')]
        execute('k8', [args.k8, args.js, *js_opts, args.alt, args.sam], orig)
        execute('rust', [args.rust, '-t', str(args.threads), *rs_opts, args.alt, args.sam], rust)
        d = first_difference(orig, rust)
        if d:
            line, x, y = d
            print(f'FAILED: difference at line {line}')
            print('k8:  ', x.rstrip()[:1000])
            print('rust:', y.rstrip()[:1000])
            raise SystemExit(1)
        if args.hla_prefix:
            lhs = {p.name[3:]: p.read_bytes() for p in Path(tmp).glob('js.*.fq')}
            rhs = {p.name[3:]: p.read_bytes() for p in Path(tmp).glob('rs.*.fq')}
            if lhs != rhs:
                print('FAILED: HLA FASTQ outputs differ')
                raise SystemExit(1)
        print('PASS: byte-for-byte identical SAM' + (' and HLA FASTQ' if args.hla_prefix else ''))


if __name__ == '__main__':
    main()
