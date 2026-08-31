#!/bin/bash

# ============================================================================
# Script: 3.Alignment.sh
# Purpose: Genome alignment and transcriptome quantification for RNA-seq reads
# Tools: STAR (Spliced Transcripts Alignment to a Reference)
# Input: Quality-filtered paired-end FASTQ files
# Output: BAM files (coordinate-sorted), transcriptome SAM files for RSEM
# ============================================================================

# Load configuration from config.sh
# Edit config.sh with your local paths before running this script
source ../config.sh

# ============================================================================
# Step 1: Generate STAR genome index
# ============================================================================
# Creates a searchable index of the human genome with splice junction information
# This step only needs to be run once per reference genome
# Runtime: ~1-2 hours depending on system resources

echo "Step 1: Creating STAR genome index..."
echo "Genome FASTA: ${GENOME_FASTA}"
echo "Annotation GTF: ${ANNOTATION_GTF}"
echo "Output index directory: ${INDEX_DIR}"

mkdir -p "${INDEX_DIR}"

STAR \
  --runThreadN "${N_THREADS}" \
  --runMode genomeGenerate \
  --genomeDir "${INDEX_DIR}" \
  --genomeFastaFiles "${GENOME_FASTA}" \
  --sjdbGTFfile "${ANNOTATION_GTF}" \
  --sjdbOverhang "${STAR_OVERHANG}"

echo "Step 1 complete: Genome index created."

# ============================================================================
# Step 2: Align reads to genome using STAR
# ============================================================================
# Align quality-filtered reads to the reference genome
# - Generates coordinate-sorted BAM files for downstream analysis
# - Generates transcriptome BAM files for RSEM quantification
# Runtime: ~30-60 minutes per sample depending on file size and resources

echo ""
echo "Step 2: Aligning reads to genome with STAR..."
echo "Index directory: ${INDEX_DIR}"
echo "Output directory: ${STAR_OUTPUT_DIR}"

mkdir -p "${STAR_OUTPUT_DIR}"

for sample_num in $(seq -f "%02g" 1 "${N_SAMPLES}"); do
  sample_name="${SAMPLE_PREFIX}${sample_num}"
  
  r1_in="${TRIMMED_DIR}/${sample_name}-R1.clean.final.fastq.gz"
  r2_in="${TRIMMED_DIR}/${sample_name}-R2.clean.final.fastq.gz"
  
  # Check if input files exist before processing
  if [ ! -f "${r1_in}" ] || [ ! -f "${r2_in}" ]; then
    echo "WARNING: Input files not found for ${sample_name}"
    echo "  Expected: ${r1_in}"
    echo "  Expected: ${r2_in}"
    continue
  fi
  
  echo "Processing: ${sample_name}"
  
  # Run STAR alignment
  STAR \
    --genomeDir "${INDEX_DIR}" \
    --runThreadN "${N_THREADS}" \
    --readFilesIn "${r1_in}" "${r2_in}" \
    --readFilesCommand zcat \
    --outFileNamePrefix "${STAR_OUTPUT_DIR}/${sample_name}." \
    --outSAMtype BAM "${STAR_BAM_SORT}" \
    --quantMode TranscriptomeSAM
    
  echo "Alignment complete for: ${sample_name}"
  
done

echo ""
echo "Step 2 complete: All samples aligned."
echo ""
echo "============================================================"
echo "Alignment pipeline complete!"
echo "============================================================"
echo "Output BAM files (coordinate-sorted): ${STAR_OUTPUT_DIR}/*.Aligned.sortedByCoord.out.bam"
echo "Output transcriptome SAM (for RSEM): ${STAR_OUTPUT_DIR}/*.Aligned.toTranscriptome.out.sam"

# ============================================================================
# END OF SCRIPT
# ============================================================================
