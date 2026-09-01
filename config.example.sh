#!/bin/bash
# ============================================================================
# Configuration Template for Bash Preprocessing Pipeline
# ============================================================================
#
# PURPOSE:
# This file demonstrates how to configure paths for the RNA-seq preprocessing
# pipeline (quality control, read trimming, alignment, quantification).
#
# INSTRUCTIONS:
# 1. Copy this file to a new file: config.sh
# 2. Edit the paths in config.sh to match your local environment
# 3. Pipeline scripts load config.sh relative to their own location
# 4. Do NOT commit config.sh to the repository
#
# ============================================================================

# ============================================================================
# PROJECT AND DATA DIRECTORIES
# ============================================================================

# Project root directory
PROJECT_ROOT="/path/to/project/root"

# Raw FASTQ input directory (contains .fastq.gz files)
DATA_INPUT_DIR="${PROJECT_ROOT}/rawdata"

# ============================================================================
# REFERENCE GENOME AND ANNOTATION
# ============================================================================

# Path to human genome FASTA file (Homo sapiens GRCh38 primary assembly)
GENOME_FASTA="/path/to/reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa"

# Path to genome annotation file (GTF format, GRCh38 release 83)
ANNOTATION_GTF="/path/to/reference/Homo_sapiens.GRCh38.83.gtf"

# Directory to store STAR genome index (will be created if doesn't exist)
INDEX_DIR="${PROJECT_ROOT}/index"

# ============================================================================
# OUTPUT DIRECTORIES (Quality Control)
# ============================================================================

# FastQC output directory
FASTQC_OUTPUT="${PROJECT_ROOT}/Results/QC/fastqc_output"

# MultiQC output directory
MULTIQC_OUTPUT="${PROJECT_ROOT}/Results/QC/multiqc_output"

# ============================================================================
# OUTPUT DIRECTORIES (Read Trimming)
# ============================================================================

# Trimmed paired-end reads output directory
TRIMMED_DIR="${PROJECT_ROOT}/Results/Preprocessing/trimmed"

# Unpaired reads output directory (from Trimmomatic)
UNPAIRED_DIR="${PROJECT_ROOT}/Results/Preprocessing/unpaired"

# ============================================================================
# OUTPUT DIRECTORIES (Alignment)
# ============================================================================

# STAR alignment output directory
STAR_OUTPUT_DIR="${PROJECT_ROOT}/Results/Alignment"

# ============================================================================
# OUTPUT DIRECTORIES (Quantification)
# ============================================================================

# RSEM quantification output directory
RSEM_OUTPUT_DIR="${PROJECT_ROOT}/Results/Quantification/RSEM"

# Dedicated RSEM reference directory and reference prefix.
# RSEM will build its own STAR index here from the same GRCh38 FASTA/GTF.
RSEM_REFERENCE_DIR="${PROJECT_ROOT}/reference/rsem"
RSEM_REFERENCE_PREFIX="${RSEM_REFERENCE_DIR}/GRCh38_rsem"

# Pipeline logs
LOG_DIR="${PROJECT_ROOT}/Results/logs"

# ============================================================================
# COMPUTATIONAL RESOURCES
# ============================================================================

# Number of CPU threads for STAR and RSEM (adjust based on your system)
N_THREADS=32

# ============================================================================
# SAMPLE INFORMATION
# ============================================================================

# Historical preprocessing used 18 technical sequencing IDs.
# These technical IDs are separate from downstream biological sample names.
N_SAMPLES=18

# Technical sample names: t01, t02, ..., t18
SAMPLE_PREFIX="t"

# ============================================================================
# TRIMMOMATIC PARAMETERS
# ============================================================================

# Quality threshold for average quality filtering (Trimmomatic AVGQUAL)
AVGQUAL=25

# Minimum read length after trimming (Trimmomatic MINLEN)
MINLEN=70

# Illumina adapter sequences for removal (TruSeq universal adapters)
ILLUMINA_ADAPTERS_R1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
ILLUMINA_ADAPTERS_R2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

# ============================================================================
# STAR PARAMETERS
# ============================================================================

# Read overhang for STAR sjdb (typically 100 for 101bp or longer reads)
STAR_OVERHANG=100

# STAR output BAM sorting method
STAR_BAM_SORT="SortedByCoordinate"

# ============================================================================
# END OF CONFIGURATION
# ============================================================================
