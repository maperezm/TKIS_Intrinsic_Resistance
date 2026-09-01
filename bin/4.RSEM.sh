#!/bin/bash

# ============================================================================
# Script: 4.RSEM.sh
# Purpose: Gene- and transcript-level RNA-seq quantification
# Tools: RSEM + STAR
#
# Input:
#   Quality-filtered paired-end FASTQ files from 2.Cleaning.sh
#
# Output:
#   RSEM expected counts, TPM and FPKM estimates
#
# RSEM runs STAR internally using an RSEM-prepared STAR reference.
# The standalone STAR workflow in 3.Alignment.sh is retained for genome
# alignment assessment and QC and is not used as RSEM input.
# ============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../config.sh"

mkdir -p "${RSEM_OUTPUT_DIR}"
mkdir -p "${RSEM_REFERENCE_DIR}"

# ============================================================================
# Step 1: Prepare RSEM + STAR reference
# ============================================================================

if [ -s "${RSEM_REFERENCE_PREFIX}.grp" ] && \
   [ -s "${RSEM_REFERENCE_DIR}/Genome" ]; then

  echo "Existing RSEM/STAR reference detected."

else

  echo "Preparing RSEM reference and STAR index..."
  echo "Genome FASTA: ${GENOME_FASTA}"
  echo "Annotation GTF: ${ANNOTATION_GTF}"

  rsem-prepare-reference \
    --gtf "${ANNOTATION_GTF}" \
    --star \
    -p "${N_THREADS}" \
    "${GENOME_FASTA}" \
    "${RSEM_REFERENCE_PREFIX}"

  echo "RSEM/STAR reference preparation complete."

fi

# ============================================================================
# Step 2: Quantification
# ============================================================================

echo ""
echo "Starting RSEM + STAR quantification..."
echo "Input directory: ${TRIMMED_DIR}"
echo "Output directory: ${RSEM_OUTPUT_DIR}"

for i in $(seq 1 "${N_SAMPLES}"); do

  sample_name=$(printf "%s%02d" "${SAMPLE_PREFIX}" "${i}")

  r1_in="${TRIMMED_DIR}/${sample_name}-R1.clean.final.fastq.gz"
  r2_in="${TRIMMED_DIR}/${sample_name}-R2.clean.final.fastq.gz"

  if [ ! -f "${r1_in}" ] || [ ! -f "${r2_in}" ]; then
    echo "WARNING: Clean FASTQ files not found for ${sample_name}"
    echo "  Expected: ${r1_in}"
    echo "  Expected: ${r2_in}"
    continue
  fi

  sample_output_dir="${RSEM_OUTPUT_DIR}/${sample_name}"
  mkdir -p "${sample_output_dir}"

  echo ""
  echo "Processing: ${sample_name}"

  rsem-calculate-expression \
    --paired-end \
    --star \
    --star-gzipped-read-file \
    -p "${N_THREADS}" \
    "${r1_in}" \
    "${r2_in}" \
    "${RSEM_REFERENCE_PREFIX}" \
    "${sample_output_dir}/${sample_name}"

  echo "Quantification complete for: ${sample_name}"

done

echo ""
echo "============================================================"
echo "RSEM + STAR quantification pipeline complete."
echo "============================================================"
echo "Output directory: ${RSEM_OUTPUT_DIR}"
echo ""
echo "Expected files per sample:"
echo "  <sample>/<sample>.genes.results"
echo "  <sample>/<sample>.isoforms.results"
echo "  <sample>/<sample>.stat/"
echo ""
echo "Downstream DGE requires biological sample identifiers and metadata."
echo "The technical-to-biological sample mapping is not included in this repository."
