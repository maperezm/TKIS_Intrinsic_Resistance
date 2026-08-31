#!/bin/bash

# ============================================================================
# Script: 4.RSEM.sh
# Purpose: Transcript-level quantification of RNA-seq reads
# Tools: RSEM (RNA-Seq by Expectation Maximization)
# Input: STAR transcriptome SAM files
# Output: Gene and transcript-level counts, FPKM values
# ============================================================================

# Load configuration from config.sh
# Edit config.sh with your local paths before running this script
source ../config.sh

# ============================================================================
# Step 1: Prepare RSEM reference (if needed)
# ============================================================================
# This step prepares the reference genome for RSEM quantification
# Only needs to be run once per reference genome
# Creates binary indices and bowtie index

echo "Note: RSEM reference preparation step"
echo "If this is your first time using this reference, you may need to run:"
echo "  rsem-prepare-reference --gtf <GTF_FILE> <GENOME_FASTA> <RSEM_REFERENCE_PREFIX>"
echo ""

# ============================================================================
# Step 2: Quantify transcripts using RSEM
# ============================================================================
# Estimates transcript and gene-level abundance from STAR transcriptome BAM
# Outputs: expected counts, FPKM, TPM, transcript counts, confidence intervals

echo "Step 1: Performing RSEM transcript quantification..."
echo "Input directory: ${STAR_OUTPUT_DIR}"
echo "Output directory: ${RSEM_OUTPUT_DIR}"

mkdir -p "${RSEM_OUTPUT_DIR}"

for sample_num in $(seq -f "%02g" 1 "${N_SAMPLES}"); do
  sample_name="${SAMPLE_PREFIX}${sample_num}"
  
  # Input: STAR transcriptome SAM file
  transcriptome_sam="${STAR_OUTPUT_DIR}/${sample_name}.Aligned.toTranscriptome.out.sam"
  
  # Check if input file exists
  if [ ! -f "${transcriptome_sam}" ]; then
    echo "WARNING: Transcriptome SAM file not found for ${sample_name}"
    echo "  Expected: ${transcriptome_sam}"
    continue
  fi
  
  echo "Processing: ${sample_name}"
  
  # Run RSEM quantification
  # Note: You must have prepared the RSEM reference index first
  # See comment in Step 1 above
  rsem-calculate-expression \
    --star \
    --paired-end \
    -p "${N_THREADS}" \
    "${transcriptome_sam}" \
    "${RSEM_REFERENCE_PREFIX}" \
    "${RSEM_OUTPUT_DIR}/${sample_name}"
  
  echo "Quantification complete for: ${sample_name}"
  
done

echo ""
echo "============================================================"
echo "RSEM transcript quantification pipeline complete!"
echo "============================================================"
echo "Output directory: ${RSEM_OUTPUT_DIR}"
echo ""
echo "Generated files per sample:"
echo "  - ${sample_name}.genes.results (gene-level quantification)"
echo "  - ${sample_name}.isoforms.results (transcript-level quantification)"
echo "  - ${sample_name}.stat/ (statistical information)"
echo ""
echo "Next step: Use these .genes.results files as input to DGE_TKIS.R"

# ============================================================================
# END OF SCRIPT
# ============================================================================
