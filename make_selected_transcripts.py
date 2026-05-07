#!/usr/bin/env python3
from __future__ import annotations

import csv
import gzip
import os
import re
import subprocess
from pathlib import Path

# The last GTF column looks like: key "value"; key "value";
# A small regex is enough to pull out gene_id, transcript_id, RPKM1, etc.
ATTR_RE = re.compile(r'(\S+) "([^"]*)"')
WORK_DIR = Path(os.environ.get('WORK_DIR', 'clean_pnas2017_fig5a'))
GTF = WORK_DIR / 'annotation/source/wgEncodeCshlLongRnaSeqGm12878CellPamTranscriptGencV7.gtf.gz'
GTF_URL = 'https://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeCshlLongRnaSeq/wgEncodeCshlLongRnaSeqGm12878CellPamTranscriptGencV7.gtf.gz'
OUT_BED = WORK_DIR / 'annotation/fig5a_selected_transcripts_before_6kb.bed'
OUT_META = WORK_DIR / 'annotation/fig5a_selected_transcripts_before_6kb.tsv'


def attrs(text):
    # Convert the GTF attribute string into a small dict.
    return dict(ATTR_RE.findall(text))


def open_text(path):
    # Read gzipped and plain text files with the same code below.
    return gzip.open(path, 'rt') if path.suffix == '.gz' else open(path)


GTF.parent.mkdir(parents=True, exist_ok=True)
# Download the GM12878 long RNA-seq transcript GTF if it is not already present.
if not GTF.exists():
    subprocess.run(['curl', '-L', '--fail', '-o', str(GTF), GTF_URL], check=True)

# First pass at transcript level: keep expressed, long, strand-annotated transcripts.
rows = []
with open_text(GTF) as inp:
    for row in csv.reader(inp, delimiter='\t'):
        if not row or row[0].startswith('#') or row[2] != 'transcript':
            continue
        chrom, start, end, strand, attr_text = row[0], int(row[3]) - 1, int(row[4]), row[6], row[8]
        a = attrs(attr_text)
        rpkm1 = float(a.get('RPKM1', 0) or 0)
        rpkm2 = float(a.get('RPKM2', 0) or 0)
        # mean(RPKM1, RPKM2) * 1000, written this way to avoid an extra division.
        score = (rpkm1 + rpkm2) * 500.0
        length = end - start
        if score >= 300 and length >= 15000 and strand in {'+', '-'}:
            rows.append((a['gene_id'], a['transcript_id'], chrom, start, end, strand, rpkm1, rpkm2, score, length))

# One gene can have many transcripts. For this gene-level profile,
# choose one representative transcript per gene before the 6 kb filter.
best = {}
for row in rows:
    gene, tx, chrom, start, end, strand, rpkm1, rpkm2, score, length = row
    old = best.get(gene)
    # Prefer higher expression, then longer transcripts; remaining ties are resolved deterministically.
    if old is None or (score, length, -start, tx) > (old[8], old[9], -old[3], old[1]):
        best[gene] = row

OUT_BED.parent.mkdir(parents=True, exist_ok=True)
selected = sorted(best.values(), key=lambda x: (x[2], x[3], x[4], x[1]))
# BED is used by the profile script; TSV keeps extra fields for checking the selection.
with open(OUT_BED, 'w', newline='') as bed, open(OUT_META, 'w', newline='') as meta:
    bed_writer = csv.writer(bed, delimiter='\t', lineterminator='\n')
    meta_writer = csv.writer(meta, delimiter='\t', lineterminator='\n')
    meta_writer.writerow(
        ['gene_id', 'transcript_id', 'chrom', 'start', 'end', 'strand', 'RPKM1', 'RPKM2', 'expression_score', 'length'])
    for gene, tx, chrom, start, end, strand, rpkm1, rpkm2, score, length in selected:
        bed_writer.writerow([chrom, start, end, tx, score, strand])
        meta_writer.writerow([gene, tx, chrom, start, end, strand, rpkm1, rpkm2, score, length])

print(len(selected), OUT_BED)
