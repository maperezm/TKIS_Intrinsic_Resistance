#function created for Diff Expression Analysis 


#---------EnhanV_plot function----
#' Generate an enhanced volcano plot
#' 
#' @param dge_obj A data frame containing differential gene expression data
#' @param title_graph Title for the volcano plot
#' @param FDR_cutoff False discovery rate cutoff for significance
#' @param point_size Size of points in the plot
#' @param titlesize Font size for the title
#' @param legendSize Font size for the legend
#' @param select_label Label for selected points (optional)
#' @param mart_version Biomart Mart object for gene information retrieval
#' 
#' @return An EnhancedVolcano plot object
#' @importFrom biomaRt getBM
#' @importFrom dplyr mutate
#' @importFrom EnhancedVolcano EnhancedVolcano
enhanced_volcano_plot <- function(dge_obj, title_graph, 
                                FDR_cutoff = 0.05,logFC_cutoff=NA,
                                  point_size = 2.5, titlesize = 22, 
                                  legendSize = 18, select_label=NA
                                ) {
  
  # Biomart query to retrieve gene information
  biomart_info <- getBM(
    attributes = c("ensembl_gene_id", 
                   "gene_biotype",
                   "external_gene_name"),
    filters    = "ensembl_gene_id",
    values     = dge_obj,
    mart       = mart_version
  )
  
  # Merge gene biotype information with dge_obj
  dge_obj <- merge(dge_obj, biomart_info, by = 1)
  
  # Define external gene names as rownames
  rownames(dge_obj) <- dge_obj$ensembl_gene_id
  
  
  # Define color based on conditions
  dge_obj <- dge_obj %>%
    mutate(color = case_when(
      FDR <= FDR_cutoff & logFC <= -logFC_cutoff & gene_biotype != "lncRNA" ~ "#6495ED",  # Blue for downregulated non-lncRNA genes
      FDR <= FDR_cutoff & logFC <= -logFC_cutoff & gene_biotype == "lncRNA" ~ "brown",  # Purple for downregulated lncRNA genes
      FDR < FDR_cutoff & logFC >= logFC_cutoff & gene_biotype != "lncRNA" ~ "#cc0000",   # Red for upregulated non-lncRNA genes
      FDR <= FDR_cutoff & logFC >= logFC_cutoff & gene_biotype == "lncRNA" ~ "purple",   # Purple for upregulated lncRNA genes
      TRUE ~ "grey50"  # Grey for other cases
    ))
  
  
  # Define key-value pairs for color labels
  keyvals_colour_y <- dge_obj$color
  keyvals_colour_y[is.na(keyvals_colour_y)] <- 'grey50'
  names(keyvals_colour_y)[keyvals_colour_y == 'grey50'] <- 'non significatives'
  names(keyvals_colour_y)[keyvals_colour_y == '#cc0000'] <- 'mRNA upregulates'
  names(keyvals_colour_y)[keyvals_colour_y == '#6495ED'] <- 'mRNA downregulates'
  names(keyvals_colour_y)[keyvals_colour_y == 'purple'] <- 'lncRNA up significatives'
  names(keyvals_colour_y)[keyvals_colour_y == 'brown'] <- 'lncRNA dow significatives'
  
  # Generate volcano plot
  volcano_p <- EnhancedVolcano(
    dge_obj,
    x = "logFC",
    y="FDR",
    lab = dge_obj$external_gene_name,
    title = title_graph,
    drawConnectors = TRUE,
    pCutoff = 0.05,
    FCcutoff = logFC_cutoff,
    pointSize = point_size,
    colAlpha = 9 / 16,
    titleLabSize = titlesize,
    caption = "",
    colCustom = keyvals_colour_y,
    #selectLab = as.character(dge_obj$external_gene_name)[which(names(keyvals_colour_y) %in% c('lncRNA up significatives', 'lncRNA dow significatives'))],  
    selectLab= select_label,
    boxedLabels = T,
    labCol = 'black',
    labFace = 'bold',
    subtitle = NULL,
    legendPosition = "bottom",
    cutoffLineWidth = 1.7,
    borderWidth = 2.3,
    legendLabSize = legendSize,
    cutoffLineCol = "black"
  )
  
  print(volcano_p)  # Print the volcano plot
}


#----------Filter RNA upregulated Genes-------
select_RNA <- function(genesDE_obj, FDR_value = NA, logFC_Value = NA) {
  
  # Input validation
  if (!is.numeric(FDR_value) || !is.numeric(logFC_Value)) {
    stop("Error: FDR_value and logFC_Value must be numeric.")
  }
  
  # Filter genes based on FDR and logFC criteria
  if (!is.na(logFC_Value)) {
    if (logFC_Value >= 1) {
      mRNA_upregulated <- genesDE_obj %>%
        filter(FDR <= FDR_value, logFC >= logFC_Value) %>%
        dplyr::select(genes, logFC, FDR)
    } else if (logFC_Value <= -1) {
      mRNA_upregulated <- genesDE_obj %>%
        filter(FDR <= FDR_value, logFC <= logFC_Value) %>%
        dplyr::select(genes, logFC, FDR)
    } else {
      stop("Error: logFC_Value should be >= 1 or <= -1.")
    }
  } else {
    mRNA_upregulated <- genesDE_obj %>%
      filter(FDR <= FDR_value) %>%
      dplyr::select(genes, logFC, FDR)
  }
  
  # Check if mRNA_upregulated is empty
  if (nrow(mRNA_upregulated) == 0) {
    message("No genes meet the specified FDR and logFC criteria.")
    return(NULL)  # Return NULL or a specific value indicating no genes found
  }
  
  # Get additional gene information from Biomart
  biomart_info <- getBM(
    attributes = c("ensembl_gene_id","external_gene_name", "entrezgene_id", "gene_biotype", "description"),
    filters = "ensembl_gene_id",
    values = mRNA_upregulated$genes,  # Use filtered gene IDs for Biomart query
    mart = mart_version  # Replace with your Biomart Mart object
  )
  
  # Merge the filtered genes with Biomart information based on ensembl_gene_id
  mRNA_upregulated_with_info <- merge(mRNA_upregulated, biomart_info, by.x = "genes", by.y = "ensembl_gene_id")
  
  return(mRNA_upregulated_with_info)
}


#------Kaplan_ meier Plot 

kaplan_meier_plot <- function( fit_obj, legend_title, data_surv)
{
  plot_1 <- ggsurvplot(fit_obj, 
                       conf.int=F, 
                       pval=T,
                       data = data_surv,
                       font.title = c(28,"bold","black"), 
                       legend.labs=c("over expressed", "under expressed"), 
                       legend.title=legend_title, 
                       legend = "right",
                       palette=c("red", "blue"), 
                       title="Kaplan-Meier Curve for Lung Cancer Survival from TCGA Dataset", 
                       break.time.by = 365,
                       xlim = c(0,2000),
                       ggtheme = theme_gray(), 
                       font.x = c(28,"bold"),
                       font.y = c(28,"bold"),
                       font.caption = c(28,"bold"), 
                       font.legend = c(27,"bold"), 
                       font.tickslab = c(28,"bold") +
                         annotate("text", label = "P-value", cex=3.3, vjust=0, hjust = 1.1, fontface=2))          
  
  return(plot_1)  
}


# Define a function to sort a dataframe by two columns and assign names to rows
sort_and_assign_names_df <- function(input_df, col1, col2, decreasing1 = TRUE, decreasing2 = TRUE) {
  sorted_df <- input_df[order(input_df[[col1]], input_df[[col2]], decreasing = c(decreasing1, decreasing2)), ]
  rownames(sorted_df) <- paste(input_df[[col1]], input_df[[col2]], sep = "_")
  return(sorted_df)
}

########
# Get additional gene information from Biomart

get_gene_info <- function(data_frame_list=NA, filter_select="ensembl_gene_id"){

biomart_info <- getBM(
  attributes = c("ensembl_gene_id",
                 "external_gene_name", 
                 "entrezgene_id", 
                 "gene_biotype", 
                 "description"),
  filters = filter_select,
  values = data_frame_list,  # Use filtered gene IDs for Biomart query
  mart = mart_version  # Replace with your Biomart Mart object
)
}



