#!/usr/bin/env python3
import gzip
import sys


def main() -> int:
    if len(sys.argv) != 3:
        sys.stderr.write('usage: dedup_length_filter_fastq_gz.py input.fastq.gz output.fastq.gz\n')
        return 2

    in_path, out_path = sys.argv[1], sys.argv[2]
    # Store only sequences here, so deduplication is sequence-based, not position-based.
    seen = set()
    raw = kept = too_short = too_long = duplicate = 0

    with gzip.open(in_path, 'rt') as inp, gzip.open(out_path, 'wt') as out:
        while True:
            # FASTQ uses four lines per read: header, sequence, plus, and quality.
            header = inp.readline()
            if not header:
                break
            seq = inp.readline()
            plus = inp.readline()
            qual = inp.readline()

            raw += 1
            seq_key = seq.strip()
            length = len(seq_key)
            # XR-seq reads should be short oligos; discard reads outside the expected range.
            if length < 20:
                too_short += 1
                continue
            if length > 50:
                too_long += 1
                continue
            # For identical sequences, keep only the first occurrence.
            if seq_key in seen:
                duplicate += 1
                continue

            seen.add(seq_key)
            out.write(header)
            out.write(seq)
            out.write(plus)
            out.write(qual)
            kept += 1

    sys.stderr.write(
        f'raw={raw} kept={kept} duplicate={duplicate} too_short={too_short} too_long={too_long}\n'
    )
    return 0


if __name__ == '__main__':
    main()
