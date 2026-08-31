#!/bin/bash

# ============================================================================
# Script: 1.Quality_control.sh
# Purpose: Quality control analysis of RNA-seq FASTQ files
# Tools: FastQC, MultiQC
# Input: Raw FASTQ files (.fastq.gz)
# Output: FastQC reports (.html, .zip), MultiQC summary
# ============================================================================

# Load configuration from config.sh
# Edit config.sh with your local paths before running this script
source ../config.sh

# Create output directories
mkdir -p "${FASTQC_OUTPUT}"
mkdir -p "${MULTIQC_OUTPUT}"

# ============================================================================
# Step 1: Run FastQC on all FASTQ files
# ============================================================================
# FastQC performs quality assessment including:
# - Per-base sequence quality
# - GC content distribution
# - Adapter content detection
# - Q30 read percentage

echo "Starting FastQC analysis on all FASTQ files..."
echo "Input directory: ${DATA_INPUT_DIR}"
echo "Output directory: ${FASTQC_OUTPUT}"

for fastq_file in "${DATA_INPUT_DIR}"/*.fastq.gz; do
  if [ -f "${fastq_file}" ]; then
    echo "Processing: $(basename "${fastq_file}")"
    fastqc -o "${FASTQC_OUTPUT}" "${fastq_file}"
  fi
done

echo "FastQC analysis complete."

# ============================================================================
# Step 2: Run MultiQC to generate aggregated summary
# ============================================================================
# MultiQC combines results from all FastQC reports into a single interactive
# HTML report for easy cross-sample comparison

echo ""
echo "Generating MultiQC summary report..."
echo "Input: ${FASTQC_OUTPUT}/*.zip"
echo "Output: ${MULTIQC_OUTPUT}/multiqc_report.html"

multiqc "${FASTQC_OUTPUT}"/*.zip -o "${MULTIQC_OUTPUT}"

echo "MultiQC report generated: ${MULTIQC_OUTPUT}/multiqc_report.html"
echo ""
echo "Quality control pipeline complete."

# ============================================================================
# END OF SCRIPT
# ============================================================================
