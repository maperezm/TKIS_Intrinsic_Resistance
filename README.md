README for Differential Expression Analysis Script
1. Overview

This script is designed for differential expression analysis of RNA sequencing data. It includes various functions for data preprocessing, filtering, visualization, and downstream analysis such as gene annotation and GO analysis. The script aims to facilitate a thorough examination of gene expression changes under different experimental conditions.

2. Requirements

Software:

R (version >= 4.0.0)
RStudio (optional but recommended for an enhanced user interface)
R Packages:

DESeq2
edgeR
Tximport
EnhancedVolcano
org.Hs.eg.db (or appropriate organism package)
clusterProfiler
BiocManager
biomaRt
dplyr
tidyverse
3. Installation

To install the necessary R packages, you can run the following commands in your R environment:

R
Copy code
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("DESeq2", "edgeR", "TXimnport", "org.Hs.eg.db", "clusterProfiler", "biomaRt"))

install.packages(c("ggplot2", "EnhancedVolcano", "VennDiagram", "dplyr", "tidyverse"))
4. Usage

4.1 Data Preparation
Ensure your RNA-seq data is formatted correctly. The expected input format is a count matrix with genes as rows and samples as columns. Metadata for the samples, including experimental conditions, should also be provided in a separate file.

4.2 Script Execution
The script includes several functions categorized into different sections. Here is an outline of the main functions and their usage
For questions or issues regarding the script, please contact:

Author: Mario Perez-Medina
Email: maperezm.medi@gmail.com   
Affiliation: Instituto Nacional de Enfermedades Respiratorias


5. License

This script is licensed under the [insert license type, e.g., MIT License]. Please see the LICENSE file for more details.

