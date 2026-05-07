#!/usr/bin/env python3
from __future__ import annotations

import csv
import os
import subprocess
from collections import defaultdict
from pathlib import Path

# Core step: use the transcript BED and repair-read BED files to compute Fig.5A profiles.
WORK_DIR = Path(os.environ.get('WORK_DIR', 'clean_pnas2017_fig5a'))
GENOME_SIZES = WORK_DIR / 'ref/hg19.chrom.sizes'
GENES_BED = WORK_DIR / 'annotation/fig5a_selected_transcripts_before_6kb.bed'
BED_DIR = WORK_DIR / 'bed'
OUT_DIR = WORK_DIR / 'profiles_simple'
TMP_DIR = WORK_DIR / 'profile_tmp_simple'
BEDTOOLS = os.environ.get('BEDTOOLS', 'bedtools')

# The three figure rows come from these three merged BED files.
SAMPLES = [
    ('BPDE_dG_GM12878_1h', BED_DIR / 'BPDE_dG_GM12878_1h.merged.bed'),
    ('CPD_NHF1_1h', BED_DIR / 'CPD_NHF1_1h.merged.bed'),
    ('PP64_NHF1_1h', BED_DIR / 'PP64_NHF1_1h.merged.bed'),
]


def read_sizes():
    # Chromosome lengths, used to keep generated windows inside valid bounds.
    with open(GENOME_SIZES) as f:
        return {row[0]: int(row[1]) for row in csv.reader(f, delimiter='\t')}


def read_genes(sizes):
    # Read representative transcripts from the previous step and apply a quick sanity filter.
    genes = []
    with open(GENES_BED) as f:
        for row in csv.reader(f, delimiter='\t'):
            chrom, start, end, name, score, strand = row[:6]
            start, end, score = int(start), int(end), float(score)
            if chrom in sizes and score >= 300 and end - start >= 15000:
                genes.append((chrom, start, end, name, score, strand))
    return sorted(genes)


def remove_neighbors(genes):
    # Strict 6 kb filter: if another selected transcript is nearby, remove both.
    # Already-removed transcripts still count as neighbors for later transcripts.
    keep = [True] * len(genes)
    by_chrom = defaultdict(list)
    for i, gene in enumerate(genes):
        by_chrom[gene[0]].append((i, gene))
    for chrom, rows in by_chrom.items():
        rows.sort(key=lambda x: x[1][1])
        for pos, (i, g) in enumerate(rows):
            g0, g1 = g[1] - 6000, g[2] + 6000
            j = pos + 1
            while j < len(rows) and rows[j][1][1] < g1:
                k, other = rows[j]
                if other[2] > g0:
                    keep[i] = False
                    keep[k] = False
                j += 1
    return [g for i, g in enumerate(genes) if keep[i]]


def windows():
    # Paper Fig.5A uses fixed 20 kb windows around TSS and TES, not scaled gene-body bins.
    out = []
    for x in range(-5000, 15000, 100):
        out.append(('TSS', x, x + 100))
    for x in range(-15000, 5000, 100):
        out.append(('TES', x, x + 100))
    return out


def anchor(gene, region):
    # TSS/TES anchors depend on transcript strand.
    chrom, start, end, name, score, strand = gene
    if region == 'TSS':
        return start if strand == '+' else end
    return end if strand == '+' else start


def complete(gene, sizes, wins):
    # Keep only transcripts whose full TSS/TES windows stay inside the chromosome.
    chrom, start, end, name, score, strand = gene
    for region, rel0, rel1 in wins:
        a = anchor(gene, region)
        g0, g1 = (a + rel0, a + rel1) if strand == '+' else (a - rel1, a - rel0)
        if g0 < 0 or g1 > sizes[chrom] or g1 <= g0:
            return False
    return True


def write_bins(genes, wins, out_path):
    # Convert transcript-relative bins into genomic BED intervals for bedtools.
    with open(out_path, 'w', newline='') as out:
        writer = csv.writer(out, delimiter='\t', lineterminator='\n')
        for chrom, start, end, name, score, strand in genes:
            for region, rel0, rel1 in wins:
                a = anchor((chrom, start, end, name, score, strand), region)
                g0, g1 = (a + rel0, a + rel1) if strand == '+' else (a - rel1, a - rel0)
                offset = 5000 if region == 'TSS' else 15000
                bin_id = (rel0 + offset) // 100
                bin_name = f'{name}|{region}|{bin_id}|{rel0}|{rel1}'
                writer.writerow([chrom, g0, g1, bin_name, bin_id, strand])


def count_lines(path):
    # Use the merged BED row count as the mapped-read denominator for this sample.
    return int(subprocess.check_output(['wc', '-l', str(path)]).split()[0])


def bedtools_count(bins, reads, out, strand_flag):
    # -S is opposite to transcript strand, counted here as TS; -s is NTS.
    # -F 0.5 requires at least half of the read to overlap the 100 bp bin.
    with open(out, 'w') as f:
        subprocess.run(
            [BEDTOOLS, 'intersect', '-a', str(bins), '-b', str(reads), '-wa', '-c', strand_flag, '-F', '0.5'], stdout=f,
            check=True)


def profile_from_counts(paths, mapped_reads, out_path):
    # Calculate per-transcript RPKM in each bin, then average across transcripts.
    sums = defaultdict(lambda: [0.0, 0, 0])
    for strand_class, path in paths:
        with open(path) as f:
            for row in csv.reader(f, delimiter='\t'):
                start, end, count = int(row[1]), int(row[2]), int(row[-1])
                name, region, bin_id, rel0, rel1 = row[3].split('|')
                key = (region, strand_class, int(bin_id), int(rel0), int(rel1), (int(rel0) + int(rel1)) / 2)
                # With 100 bp bins, this is count per 0.1 kb per million mapped reads.
                rpkm = count * 1_000_000_000.0 / (mapped_reads * (end - start))
                sums[key][0] += rpkm
                sums[key][1] += count
                sums[key][2] += 1
    with open(out_path, 'w', newline='') as out:
        writer = csv.writer(out, delimiter='\t', lineterminator='\n')
        writer.writerow(
            ['region', 'strand_class', 'bin_index', 'rel_start', 'rel_end', 'rel_mid', 'mean_RPKM', 'sum_counts',
             'n_genes'])
        for key in sorted(sums):
            rpkm_sum, count_sum, n = sums[key]
            writer.writerow([*key, rpkm_sum / n, count_sum, n])


OUT_DIR.mkdir(parents=True, exist_ok=True)
TMP_DIR.mkdir(parents=True, exist_ok=True)
sizes = read_sizes()
wins = windows()
# Final transcript set for plotting: representative transcripts + 6 kb filter + complete windows.
genes = [g for g in remove_neighbors(read_genes(sizes)) if complete(g, sizes, wins)]
with open(OUT_DIR / 'selected_transcripts_after_6kb.tsv', 'w', newline='') as out:
    writer = csv.writer(out, delimiter='\t', lineterminator='\n')
    writer.writerow(['chrom', 'start', 'end', 'transcript_id', 'expression_score', 'strand', 'length'])
    for chrom, start, end, name, score, strand in genes:
        writer.writerow([chrom, start, end, name, score, strand, end - start])

bins = TMP_DIR / 'fig5a_bins.bed'
write_bins(genes, wins, bins)
for label, reads in SAMPLES:
    # For each damage type, create TS/NTS count tables and convert them to one profile.
    mapped = count_lines(reads)
    ts = TMP_DIR / f'{label}.TS.counts.bed'
    nts = TMP_DIR / f'{label}.NTS.counts.bed'
    bedtools_count(bins, reads, ts, '-S')
    bedtools_count(bins, reads, nts, '-s')
    profile_from_counts([('TS', ts), ('NTS', nts)], mapped, OUT_DIR / f'{label}.profile.tsv')
    ts.unlink()
    nts.unlink()
    print(label, mapped)
bins.unlink()
print(len(genes), OUT_DIR / 'selected_transcripts_after_6kb.tsv')
