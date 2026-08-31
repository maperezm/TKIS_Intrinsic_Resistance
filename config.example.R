# ============================================================================
# Configuration Template for TKIS Intrinsic Resistance Analysis
# ============================================================================
# 
# PURPOSE:
# This file demonstrates how to configure local and external file paths
# for the RNA-seq differential expression analysis pipeline.
#
# INSTRUCTIONS:
# 1. Copy this file to a new file: config.R
# 2. Edit the paths in config.R to match your local environment
# 3. Source config.R at the beginning of DGE_TKIS.R:
#    source("config.R")
# 4. Do NOT commit config.R to the repository (it is in .gitignore)
#
# ============================================================================

# Project root directory (relative path from repository root)
# Change this if your working directory differs
PROJECT_ROOT <- here::here()

# ============================================================================
# DATA PATHS
# ============================================================================

# Path to RSEM quantification results
# RSEM output should be organized as: RSEM_DIR/sample_name/sample_name.genes.results
RSEM_DIR <- "/path/to/external/RSEM/results"

# Example: If your RSEM files are stored in a Google Drive:
# RSEM_DIR <- "~/My\ Drive/RNAseq/Analisis_Expresion/RSEM"

# Path to sample metadata file
SAMPLE_METADATA_FILE <- file.path(PROJECT_ROOT, "Data", "sample_Persistent.txt")

# ============================================================================
# REFERENCE DATA PATHS (used in bash scripts)
# ============================================================================

# Path to human genome FASTA file (GRCh38)
GENOME_FASTA <- "/path/to/reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa"

# Path to genome annotation GTF file (GRCh38)
ANNOTATION_GTF <- "/path/to/reference/Homo_sapiens.GRCh38.83.gtf"

# Path to STAR index directory (will be created if doesn't exist)
STAR_INDEX_DIR <- "/path/to/STAR/index"

# ============================================================================
# INPUT DATA PATHS (bash preprocessing pipeline)
# ============================================================================

# Raw FASTQ file directory
FASTQ_INPUT_DIR <- "/path/to/rawdata/FASTQ"

# Output directory for quality control (FastQC)
FASTQC_OUTPUT_DIR <- file.path(PROJECT_ROOT, "Results", "QC", "fastqc_output")

# Output directory for MultiQC
MULTIQC_OUTPUT_DIR <- file.path(PROJECT_ROOT, "Results", "QC", "multiqc_output")

# Directory for trimmed reads
TRIMMED_OUTPUT_DIR <- file.path(PROJECT_ROOT, "Results", "Preprocessing", "trimmed")

# Directory for unpaired reads
UNPAIRED_OUTPUT_DIR <- file.path(PROJECT_ROOT, "Results", "Preprocessing", "unpaired")

# Directory for STAR alignment outputs
STAR_OUTPUT_DIR <- file.path(PROJECT_ROOT, "Results", "Alignment", "STAR")

# Directory for RSEM outputs
RSEM_OUTPUT_DIR <- file.path(PROJECT_ROOT, "Results", "Quantification", "RSEM")

# ============================================================================
# R ANALYSIS OUTPUT PATHS
# ============================================================================

# Directory for R analysis results (figures, tables)
R_RESULTS_DIR <- file.path(PROJECT_ROOT, "Results", "DifferentialExpression")

# ============================================================================
# ANALYSIS PARAMETERS
# ============================================================================

# CPM (Counts Per Million) filtering threshold
CPM_THRESHOLD <- 1
MIN_SAMPLES_CPM <- 3  # Minimum number of samples exceeding CPM threshold

# TMM normalization method (use "TMM" for Trimmed Mean of M-values)
NORMALIZATION_METHOD <- "TMM"

# Fold-change threshold (log2 scale)
LFC_THRESHOLD <- 1

# False Discovery Rate (FDR) threshold
FDR_THRESHOLD <- 0.05

# Multiple testing correction method
CORRECTION_METHOD <- "BH"  # Benjamini-Hochberg

# ============================================================================
# BIOMART CONFIGURATION
# ============================================================================

# Biomart database version and dataset for human genes
BIOMART_SPECIES <- "hsapiens"
BIOMART_DATASET <- "hsapiens_gene_ensembl"

# ============================================================================
# SAMPLE SPECIFICATIONS
# ============================================================================

# Sample information can be modified here if needed
# Otherwise, samples are loaded from SAMPLE_METADATA_FILE

# Number of threads for parallel processing (adjust based on your system)
N_THREADS <- 24

# ============================================================================
# OUTPUT FIGURE PARAMETERS
# ============================================================================

# DPI for PNG figures
PNG_DPI <- 200

# Default figure dimensions (width, height in inches)
FIG_WIDTH_DEFAULT <- 12
FIG_HEIGHT_DEFAULT <- 8

# ============================================================================
# END OF CONFIGURATION
# ============================================================================
