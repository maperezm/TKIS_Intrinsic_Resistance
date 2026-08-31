#!/bin/bash

# ============================================================================
# Script: 2.Cleaning.sh
# Purpose: Read quality trimming and adapter removal for RNA-seq FASTQ files
# Tools: Trimmomatic, cutadapt
# Input: Raw paired-end FASTQ files (.fastq.gz)
# Output: Quality-trimmed paired and unpaired FASTQ files
# ============================================================================

# Load configuration from config.sh
# Edit config.sh with your local paths before running this script
source ../config.sh

# Create output directories
mkdir -p "${TRIMMED_DIR}"
mkdir -p "${UNPAIRED_DIR}"

# ============================================================================
# Step 1: Quality filtering with Trimmomatic (AVGQUAL)
# ============================================================================
# Filter reads based on average quality score
# - AVGQUAL:25 keeps reads with average quality >= 25
# - Removes reads falling below quality threshold
# - Output: paired and unpaired reads

echo "Step 1: Quality filtering with Trimmomatic (AVGQUAL:${AVGQUAL})..."

for sample_num in $(seq -f "%02g" 1 "${N_SAMPLES}"); do
  sample_name="${SAMPLE_PREFIX}${sample_num}"
  
  r1_input="${DATA_INPUT_DIR}/${sample_name}-R1.fastq.gz"
  r2_input="${DATA_INPUT_DIR}/${sample_name}-R2.fastq.gz"
  
  r1_paired="${TRIMMED_DIR}/${sample_name}-R1.paired.fastq.gz"
  r2_paired="${TRIMMED_DIR}/${sample_name}-R2.paired.fastq.gz"
  r1_unpaired="${UNPAIRED_DIR}/${sample_name}-R1.unpaired.fastq.gz"
  r2_unpaired="${UNPAIRED_DIR}/${sample_name}-R2.unpaired.fastq.gz"
  
  echo "Processing: ${sample_name}"
  
  # Run Trimmomatic with AVGQUAL filter
  trimmomatic PE \
    "${r1_input}" "${r2_input}" \
    "${r1_paired}" "${r1_unpaired}" \
    "${r2_paired}" "${r2_unpaired}" \
    AVGQUAL:${AVGQUAL} \
    2>&1 | tee -a "${DATA_INPUT_DIR}/quality_trimming_log.txt"
    
done

echo "Step 1 complete: Quality filtering finished."

# ============================================================================
# Step 2: Illumina adapter removal with cutadapt
# ============================================================================
# Remove 3' adapter sequences from quality-filtered reads
# - Removes TruSeq adapter sequences
# - Processes paired-end reads with -p flag to maintain read pairing

echo ""
echo "Step 2: Illumina adapter removal with cutadapt..."

for sample_num in $(seq -f "%02g" 1 "${N_SAMPLES}"); do
  sample_name="${SAMPLE_PREFIX}${sample_num}"
  
  r1_in="${TRIMMED_DIR}/${sample_name}-R1.paired.fastq.gz"
  r2_in="${TRIMMED_DIR}/${sample_name}-R2.paired.fastq.gz"
  
  r1_clean="${TRIMMED_DIR}/${sample_name}-R1.clean.paired.fastq.gz"
  r2_clean="${TRIMMED_DIR}/${sample_name}-R2.clean.paired.fastq.gz"
  
  echo "Processing: ${sample_name} - Adapter removal"
  
  # Run cutadapt for adapter trimming
  cutadapt \
    -a "${ILLUMINA_ADAPTERS_R1}" \
    -A "${ILLUMINA_ADAPTERS_R2}" \
    -o "${r1_clean}" \
    -p "${r2_clean}" \
    "${r1_in}" "${r2_in}" \
    2>&1 | tee -a "${DATA_INPUT_DIR}/adapter_trimming_log.txt"
    
done

echo "Step 2 complete: Adapter removal finished."

# ============================================================================
# Step 3: Minimum length filtering with Trimmomatic (MINLEN)
# ============================================================================
# Remove reads shorter than MINLEN threshold after all trimming steps
# - MINLEN:70 removes reads < 70 bp (typical for RNA-seq with 101+ bp reads)
# - Prevents alignment issues with very short reads

echo ""
echo "Step 3: Minimum length filtering with Trimmomatic (MINLEN:${MINLEN})..."

for sample_num in $(seq -f "%02g" 1 "${N_SAMPLES}"); do
  sample_name="${SAMPLE_PREFIX}${sample_num}"
  
  r1_in="${TRIMMED_DIR}/${sample_name}-R1.clean.paired.fastq.gz"
  r2_in="${TRIMMED_DIR}/${sample_name}-R2.clean.paired.fastq.gz"
  
  r1_final="${TRIMMED_DIR}/${sample_name}-R1.clean.final.fastq.gz"
  r2_final="${TRIMMED_DIR}/${sample_name}-R2.clean.final.fastq.gz"
  r1_unpaired_final="${UNPAIRED_DIR}/${sample_name}-R1.clean.unpaired.fastq.gz"
  r2_unpaired_final="${UNPAIRED_DIR}/${sample_name}-R2.clean.unpaired.fastq.gz"
  
  echo "Processing: ${sample_name} - Minimum length filtering"
  
  # Run Trimmomatic with MINLEN filter
  trimmomatic PE \
    "${r1_in}" "${r2_in}" \
    "${r1_final}" "${r1_unpaired_final}" \
    "${r2_final}" "${r2_unpaired_final}" \
    MINLEN:${MINLEN} \
    2>&1 | tee -a "${DATA_INPUT_DIR}/minlen_filtering_log.txt"
    
done

echo "Step 3 complete: Minimum length filtering finished."
echo ""
echo "============================================================"
echo "Read cleaning pipeline complete!"
echo "============================================================"
echo "Quality-filtered reads are in: ${TRIMMED_DIR}"
echo "Unpaired reads are in: ${UNPAIRED_DIR}"
echo "Log files:"
echo "  - ${DATA_INPUT_DIR}/quality_trimming_log.txt"
echo "  - ${DATA_INPUT_DIR}/adapter_trimming_log.txt"
echo "  - ${DATA_INPUT_DIR}/minlen_filtering_log.txt"

# ============================================================================
# END OF SCRIPT
# ============================================================================
