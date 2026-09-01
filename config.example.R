# ============================================================================
# Configuration template for the differential-expression analysis
# ============================================================================
#
# Copy this file to:
#
#   config.R
#
# Edit local paths in config.R. The real config.R is ignored by Git.
# ============================================================================

# Project root.
# When DGE_TKIS.R is run with:
#
#   Rscript bin/DGE_TKIS.R
#
# PROJECT_ROOT_DEFAULT is resolved automatically from the script location.
PROJECT_ROOT <- if (exists("PROJECT_ROOT_DEFAULT")) {
  PROJECT_ROOT_DEFAULT
} else {
  normalizePath(".", mustWork = TRUE)
}

# ============================================================================
# INPUT DATA
# ============================================================================

# Directory containing RSEM outputs organized as:
#
# RSEM_DIR/
#   HCC827_OSI_1/
#     HCC827_OSI_1.genes.results
#   ...
#
# These are biological sample identifiers used in the downstream analysis.
# Historical technical sequencing IDs (t01, t02, ...) were used during
# preprocessing. The technical-to-biological mapping is not included here.
RSEM_DIR <- "/path/to/RSEM/results"

SAMPLE_METADATA_FILE <- file.path(
  PROJECT_ROOT,
  "Data",
  "sample_Persistent.txt"
)

# ============================================================================
# OUTPUT
# ============================================================================

R_RESULTS_DIR <- file.path(PROJECT_ROOT, "Results")

# ============================================================================
# ANALYSIS PARAMETERS
# ============================================================================

# Preserve the original filtering rule:
# log2-CPM > 1 in at least 3 samples.
LOG_CPM_THRESHOLD <- 1
MIN_SAMPLES_LOG_CPM <- 3

# edgeR normalization
NORMALIZATION_METHOD <- "TMM"

# Differential-expression thresholds
LFC_THRESHOLD <- 1
FDR_THRESHOLD <- 0.05
CORRECTION_METHOD <- "BH"

# ============================================================================
# BIOMART
# ============================================================================

BIOMART_MART <- "ENSEMBL_MART_ENSEMBL"
BIOMART_DATASET <- "hsapiens_gene_ensembl"
