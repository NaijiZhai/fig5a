#!/usr/bin/env bash
set -euo pipefail

# Set paths first. For a new run, WORK_DIR is usually the only thing to change.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
WORK_DIR="${WORK_DIR:-${ROOT_DIR}/clean_pnas2017_fig5a}"
ENV_DIR="${ENV_DIR:-${CONDA_PREFIX:-/Users/zhainaiji/miniforge3/envs/pnas2017_fig5a}}"
THREADS="${THREADS:-8}"

# Prefer tools from the conda env, so system versions do not get mixed in.
export PATH="${ENV_DIR}/bin:${PATH}"

# Keep intermediate files in separate folders, so failed steps are easier to trace.
REF_DIR="${WORK_DIR}/ref"
SRA_DIR="${WORK_DIR}/sra_simple"
FASTQ_DIR="${WORK_DIR}/fastq_simple"
TRIM_DIR="${WORK_DIR}/trimmed_simple"
BAM_DIR="${WORK_DIR}/bam_simple"
BED_DIR="${WORK_DIR}/bed"
LOG_DIR="${WORK_DIR}/logs_simple"
mkdir -p "${REF_DIR}" "${SRA_DIR}" "${FASTQ_DIR}" "${TRIM_DIR}" "${BAM_DIR}" "${BED_DIR}" "${LOG_DIR}"

ADAPTER="${REF_DIR}/xrseq_3prime_adapter.fa"
HG19_GZ="${REF_DIR}/hg19.fa.gz"
HG19_FA="${REF_DIR}/hg19.fa"
HG19_FEMALE_FA="${REF_DIR}/hg19_female.fa"
HG19_INDEX="${REF_DIR}/hg19_female"

# Adapter sequence found from tail k-mers in the BPDE raw FASTQ files.
cat > "${ADAPTER}" <<'EOF'
>XRseq_3prime_adapter_core
TGGAATTCTCGGGTGCCAAGGAACTCCAGT
EOF

# Chrom sizes are used later to keep TSS/TES windows inside chromosome bounds.
if [ ! -s "${REF_DIR}/hg19.chrom.sizes" ]; then
    curl -L --fail -o "${REF_DIR}/hg19.chrom.sizes" \
        https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes
fi

# Building the Bowtie index is slow, so skip it when the index already exists.
if [ ! -s "${HG19_INDEX}.1.ebwt" ]; then
    curl -L --fail -o "${HG19_GZ}" \
        https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
    gzip -dc "${HG19_GZ}" > "${HG19_FA}"
    samtools faidx "${HG19_FA}"
    : > "${HG19_FEMALE_FA}"
    # The paper used a human female genome; here I approximate it as hg19 without chrY.
    cut -f1 "${HG19_FA}.fai" | grep -v '^chrY$' | while read -r chrom; do
        samtools faidx "${HG19_FA}" "${chrom}" >> "${HG19_FEMALE_FA}"
    done
    samtools faidx "${HG19_FEMALE_FA}"
    bowtie-build "${HG19_FEMALE_FA}" "${HG19_INDEX}"
    # The raw FASTA is large, so remove it after the index is built.
    rm -f "${HG19_GZ}" "${HG19_FA}" "${HG19_FA}.fai" "${HG19_FEMALE_FA}" "${HG19_FEMALE_FA}.fai"
fi

process_run() {
    local run="$1"
    local bed="${BED_DIR}/${run}.bed"
    local raw="${FASTQ_DIR}/${run}.fastq"
    local raw_gz="${raw}.gz"
    local trimmed="${TRIM_DIR}/${run}.trimmed.fastq.gz"
    local filtered="${TRIM_DIR}/${run}.filtered.dedup.fastq.gz"
    local bam="${BAM_DIR}/${run}.sorted.bam"

    # If this run already has a BED file, skip it for easy resume.
    if [ -s "${bed}" ]; then
        return
    fi

    # Download the SRA file, then convert it to single-end FASTQ.
    prefetch "${run}" --output-directory "${SRA_DIR}"
    fasterq-dump "${SRA_DIR}/${run}/${run}.sra" --outdir "${FASTQ_DIR}" --threads "${THREADS}" --split-files
    if [ -s "${FASTQ_DIR}/${run}_1.fastq" ]; then
        mv "${FASTQ_DIR}/${run}_1.fastq" "${raw}"
    fi
    gzip -f "${raw}"

    # Trim adapters before length filtering, otherwise adapter tails can distort read length.
    trimmomatic SE -threads "${THREADS}" -phred33 \
        "${raw_gz}" "${trimmed}" \
        "ILLUMINACLIP:${ADAPTER}:2:30:10" \
        MINLEN:1

    # Simple filter: keep 20-50 nt reads and deduplicate by read sequence.
    "${ENV_DIR}/bin/python" "${SCRIPT_DIR}/dedup_length_filter_fastq_gz.py" \
        "${trimmed}" "${filtered}"

    # Bowtie parameters mostly follow the paper; pipe into samtools to avoid extra large files.
    gzip -dc "${filtered}" \
        | bowtie -q --nomaqround --phred33-quals -m 4 -n 2 -e 70 -l 20 \
            --best -p "${THREADS}" --seed 123 -S "${HG19_INDEX}" - \
            2> "${LOG_DIR}/${run}.bowtie.log" \
        | samtools view -bS -F 4 - \
        | samtools sort -@ "${THREADS}" -o "${bam}" -

    # Fig.5A only needs BED later, so the BAM can be removed after conversion.
    bedtools bamtobed -i "${bam}" | LC_ALL=C sort -k1,1 -k2,2n > "${bed}"
    rm -rf "${SRA_DIR:?}/${run}" "${raw_gz}" "${trimmed}" "${filtered}" "${bam}"
}

merge_runs() {
    local label="$1"
    shift
    local merged="${BED_DIR}/${label}.merged.bed"
    # If the replicate-merged BED already exists, skip this group.
    if [ -s "${merged}" ]; then
        return
    fi
    # Make sure each SRR has its own BED before merging.
    for run in "$@"; do
        process_run "${run}"
    done
    # Paper Fig.5A uses merged replicates; here they are concatenated and sorted.
    : > "${merged}.tmp"
    for run in "$@"; do
        cat "${BED_DIR}/${run}.bed" >> "${merged}.tmp"
    done
    LC_ALL=C sort -k1,1 -k2,2n "${merged}.tmp" > "${merged}"
    rm -f "${merged}.tmp"
    for run in "$@"; do
        rm -f "${BED_DIR}/${run}.bed"
    done
}

# Three damage types, matching the three rows in paper Fig.5A.
merge_runs BPDE_dG_GM12878_1h SRR5444677 SRR5444678
merge_runs CPD_NHF1_1h SRR5444679 SRR5444680 SRR5444681 SRR5444682
merge_runs PP64_NHF1_1h SRR5444683 SRR5444684 SRR5444685 SRR5444686

# Print read counts at the end; these counts are used as RPKM denominators.
wc -l "${BED_DIR}/BPDE_dG_GM12878_1h.merged.bed"
wc -l "${BED_DIR}/CPD_NHF1_1h.merged.bed"
wc -l "${BED_DIR}/PP64_NHF1_1h.merged.bed"
