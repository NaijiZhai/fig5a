# Fig.5A workflow

This folder is the simplified pipeline for reproducing a Fig.5A-like TSS/TES repair profile.

Scripts:

- `make_reads_bed.sh`: SRR -> FASTQ -> trim adapter -> length/dedup filter -> bowtie -> merged BED
- `dedup_length_filter_fastq_gz.py`: small FASTQ length filter and sequence-level dedup helper
- `make_selected_transcripts.py`: GM12878 RNA-seq GTF -> selected representative transcripts
- `make_profiles.py`: selected transcripts + repair BED -> TS/NTS RPKM profiles
- `plot_fig5a.py`: profile TSV -> figure PNG/PDF

Main choices:

- `expression_score = mean(RPKM1, RPKM2) * 1000`
- expression cutoff `>= 300`
- transcript length `>= 15000`
- one representative transcript per gene
- strict neighbor exclusion `+/- 6000 bp`, any strand
- TSS window `-5 kb` to `+15 kb`
- TES window `-15 kb` to `+5 kb`
- bin size `100 bp`
- TS: `bedtools intersect -S -F 0.5`
- NTS: `bedtools intersect -s -F 0.5`

Create conda environment first:

```bash
conda create -n pnas2017_fig5a -c bioconda -c conda-forge \
  python=3.10 sra-tools trimmomatic bowtie samtools bedtools matplotlib -y
```

Run from repo root:

```bash
conda activate pnas2017_fig5a

export WORK_DIR="$PWD/fig5a"
export ENV_DIR="$CONDA_PREFIX"
export THREADS=8

bash make_reads_bed.sh
python make_selected_transcripts.py
BEDTOOLS="${ENV_DIR}/bin/bedtools" python make_profiles.py
python plot_fig5a.py
```

Main outputs:

```text
${WORK_DIR}/bed/BPDE_dG_GM12878_1h.merged.bed
${WORK_DIR}/bed/CPD_NHF1_1h.merged.bed
${WORK_DIR}/bed/PP64_NHF1_1h.merged.bed
${WORK_DIR}/profiles_simple/selected_transcripts_after_6kb.tsv
${WORK_DIR}/profiles_simple/fig5a_raw.png
${WORK_DIR}/profiles_simple/fig5a_paper_like_x1.6.png
```
