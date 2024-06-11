##Script designed to perform DE analysis from persistent cells 
# Check if a package is installed and install it if not
check_and_install <- function(package){
  if (!require(package, character.only = TRUE)){
    install.packages(package)
    library(package, character.only = TRUE)
  }
}

# List of packages
packages <- c("edgeR"       ,  "tximport", 
              "ggplot2"     ,  "gridExtra", 
              "tidyverse"   ,  "RVenn", 
              "biomaRt"     ,  "EnhancedVolcano", 
              "pheatmap"    ,  "PCAtools",
              "clusterProfiler", "RColorBrewer", 
              "org.Hs.eg.db", "here","enrichplot")

# Check and install packages
lapply(packages, check_and_install)

# Set working directory
here <- here("~/Dropbox/Trabajo/Repositorios/Intrinsicr_Resistance/TKIS/") #principal path directory where you put dowload the files

setwd(here)

# Source functions
source("bin/functions.R")

# Use Mart
mart_version <-  useMart( biomart='ENSEMBL_MART_ENSEMBL', dataset='hsapiens_gene_ensembl')

# Set RSEM_directory
dir <- ("~/maperezm.medi@gmail.com - Google Drive/My Drive/Doctorado_Driv/RNAseq/Analisis_Expresion/RSEM") #path you put RSEM files and TXT with description files

# Read samples
samples <- read.table(file.path(here, "Data/sample_Persistent.txt"), he =T)

# Extract conditions and cell line information
samples <- mutate(samples, 
   Condition = str_extract(samples$sample, "CONTROL|ERLO|OSI|DMSO"), 
   Cell_line = str_extract(samples$sample, "HCC827|HCC4006|H1975"))

# Construct file paths
files_genes <- file.path(dir, samples$sample, paste0(samples$sample, ".genes.results"))

# Check if files exist
missing_files <- files_genes[!file.exists(files_genes)]
if (length(missing_files) > 0) {
  cat("Los siguientes archivos no están presentes:\n")
  cat(paste("- ", missing_files, "\n", sep = ""), sep = "")
  }

# Load files using tximport
txi.rsem_genes <- tximport(files_genes, 
                           type         = "rsem", 
                           geneIdCol    = "gene_id",
                           txIdCol      = "transcript_id", 
                           countsCol    = "expected_count",
                           lengthCol    = "effective_length", 
                           abundanceCol = "FPKM")

# Convert counts to integer and store as dataframe
counts <- txi.rsem_genes$counts
storage.mode(counts) <- "integer"
counts <- as.data.frame(counts)
colnames(counts) <- samples$sample
rows_with_na <- rowSums(is.na(counts)) > 0
counts <- counts[!rows_with_na, ]


# Filter genes with low abundance
genes_pass_filter <- rowSums(cpm(counts, log = T) > 1) > 2

##Visualize how many genes passed the filtering criteria
table(genes_pass_filter)

##keep the genes which passed filtering criteria in the original matrix
counts <- counts[genes_pass_filter, ]

##Create a factor object to indicate the experimental conditions
groups <- factor(sub("..$", "", names(counts))) ##"..$" indicates substitute the two last characters with "" (nothing)
table(groups)

##Store the counts and the groups in an edgeR object
edgeRlist <- DGEList(counts = counts, 
      group = groups, 
      genes = rownames(counts))


##Calculate the normalization factors using the TMM method
edgeRlist <- calcNormFactors(edgeRlist, method = "TMM")

##Visualize the normalization factors
edgeRlist$samples

row.names(samples) <- samples$sample


##Plot the results using relative expression in each sample
pdf("Results/mdplot_all.pdf", width = 10, height = 10, pointsize = 12)
par(mfrow = c(2, 2)) # 2 rows, 2 columns
for (i in 1:ncol(edgeRlist)) {
  plotMD(cpm(edgeRlist, log = TRUE), column = i, 
         xlab = "Average log-expression", 
         ylab = "Expression log-ratio (this sample vs others)", 
         main = colnames(edgeRlist)[i])
  abline(h = 0, col = "red", lty = 2, lwd = 2)
}
dev.off()

##Inspect replicates by performing a PCA analysis
pca <- pca(cpm(edgeRlist$counts, log = T,
               normalized.lib.sizes = T), 
           scale= T, metadata= samples)

pdf("Results/PCAplot.pdf", width = 14, height = 10.2)
biplot(pca, 
       lab = colnames(edgeRlist$counts), 
       pointSize = 5, 
       title = "PCA", 
       labSize = 5,
       titleLabSize = 20, colby = "Condition",
       legendPosition = "bottom", axisLabSize = 15,
       legendLabSize = 15, colkey = c("#168c0d",
                                      "#daa920",
                                      "#d947a9",
                                      "lightblue",
                                      "red"))   
dev.off()
       
#----------------------Differential expression analysis--------------------
##Make the design matrix
design <- model.matrix(~0+edgeRlist$samples$group)


##Asign colnames of the design matrix by using the levels of the experimental groups
colnames(design) <- levels(edgeRlist$samples$group)

##Estimate data dispersion
edgeRlist <- estimateDisp(edgeRlist, design = design, robust = T)

##Visualize the dispersion levels
png("Results/data_dispersion.png", height = 700, width = 800)
plotBCV(edgeRlist)
dev.off()

##Get the contrast matrix compared residual cells vs Control cells 
contrast <- makeContrasts("HCC827_ERLO"  = "HCC827_ERLO  - HCC827_CONTROL" ,
                          "HCC827_DMSO"  = "HCC827_DMSO  - HCC827_CONTROL" ,
                          "HCC4006_ERLO" = "HCC4006_ERLO - HCC4006_CONTROL", 
                          "HCC4006_DMSO" = "HCC4006_DMSO - HCC4006_CONTROL",
                          "H1975_OSI"    = "H1975_OSI    - H1975_CONTROL"  , 
                          "HCC827_OSI"   = "HCC827_OSI   - HCC827_CONTROL" ,
                          "HCC827_DMSO"  = "HCC827_DMSO  - HCC827_CONTROL" ,
                          "HCC4006_OSI"  = "HCC4006_OSI  - HCC4006_CONTROL", 
                          "HCC4006_DMSO" = "HCC4006_DMSO - HCC4006_CONTROL", 
                          levels = edgeRlist$design)

##Adjust data to a negative bi-nomial generalized linear model
fit <- glmQLFit(edgeRlist, 
                design = design, robust = T,  
                dispersion = edgeRlist$trended.dispersion)


##Test the null hypothesis in which lfc of genes in the residual cells are zero respect to Control cells # nolint
geneQLF_HCC4006_ERLOvsControl <- glmQLFTest(fit, coef = ncol(fit$design), contrast = contrast[, "HCC4006_ERLO"])
geneQLF_HCC827_ERLOvsControl  <- glmQLFTest(fit, coef = ncol(fit$design), contrast = contrast[, "HCC827_ERLO"])
geneQLF_H1975_OSIvsControl    <- glmQLFTest(fit, coef = ncol(fit$design), contrast = contrast[, "H1975_OSI"])
geneQLF_HCC827_OSIvsControl   <- glmQLFTest(fit, coef = ncol(fit$design), contrast = contrast[, "HCC827_OSI"])
geneQLF_HCC4006_OSIvsControl  <- glmQLFTest(fit, coef = ncol(fit$design), contrast = contrast[, "HCC4006_OSI"])
geneQLF_HCC4006_DMSOvsControl <- glmQLFTest(fit, coef = ncol(fit$design), contrast = contrast[, "HCC4006_DMSO"])
geneQLF_HCC827_DMSOvsControl  <- glmQLFTest(fit, coef = ncol(fit$design), contrast = contrast[, "HCC827_DMSO"])
geneQLF_HCC827_DMSOvsControl  <- glmQLFTest(fit, coef = ncol(fit$design), contrast = contrast[, "HCC827_DMSO"])
geneQLF_HCC4006_DMSOvsControl <- glmQLFTest(fit, coef = ncol(fit$design), contrast = contrast[, "HCC4006_DMSO"])


#Obtain DGE with lFC 1 and FDR = 0.05 
is.de.gene_HCC4006_ERLOvsControl <- decideTests(geneQLF_HCC4006_ERLOvsControl ,  adjust.method = "BH", lfc = 1, p.value = 0.05)
is.de.gene_HCC827_ERLOvsControl  <- decideTests(geneQLF_HCC827_ERLOvsControl  ,  adjust.method = "BH", lfc = 1, p.value = 0.05)
is.de.gene_H1975_OSIvsControl    <- decideTests(geneQLF_H1975_OSIvsControl    ,  adjust.method = "BH", lfc = 1, p.value = 0.05)
is.de.gene_HCC827_OSIvsControl   <- decideTests(geneQLF_HCC827_OSIvsControl   ,  adjust.method = "BH", lfc = 1, p.value = 0.05)
is.de.gene_HCC4006_OSIvsControl  <- decideTests(geneQLF_HCC4006_OSIvsControl  ,  adjust.method = "BH", lfc = 1, p.value = 0.05)
is.de.gene_HCC4006_DMSOvsControl <- decideTests(geneQLF_HCC4006_DMSOvsControl ,  adjust.method = "BH", lfc = 1, p.value = 0.05)
is.de.gene_HCC827_DMSOvsControl  <- decideTests(geneQLF_HCC827_DMSOvsControl  ,  adjust.method = "BH", lfc = 1, p.value = 0.05)
is.de.gene_HCC827_DMSOvsControl  <- decideTests(geneQLF_HCC827_DMSOvsControl  ,  adjust.method = "BH", lfc = 1, p.value = 0.05)
is.de.gene_HCC4006_DMSOvsControl <- decideTests(geneQLF_HCC4006_DMSOvsControl ,  adjust.method = "BH", lfc = 1, p.value = 0.05)


##Visualize how many genes rejected the null hypothesis with a FDR threshold of 0.05
summary (is.de.gene_HCC4006_ERLOvsControl )
summary (is.de.gene_HCC827_ERLOvsControl  )
summary (is.de.gene_H1975_OSIvsControl    )
summary (is.de.gene_HCC827_OSIvsControl   )
summary (is.de.gene_HCC4006_OSIvsControl  )
summary (is.de.gene_HCC4006_DMSOvsControl )
summary (is.de.gene_HCC827_DMSOvsControl  )


##Store the results in a data frame
geneDE_HCC4006_ERLOvsControl  <- as.data.frame(topTags(geneQLF_HCC4006_ERLOvsControl , n = nrow(geneQLF_HCC4006_ERLOvsControl )))
geneDE_HCC827_ERLOvsControl   <- as.data.frame(topTags(geneQLF_HCC827_ERLOvsControl  , n = nrow(geneQLF_HCC827_ERLOvsControl  )))
geneDE_H1975_OSIvsControl     <- as.data.frame(topTags(geneQLF_H1975_OSIvsControl    , n = nrow(geneQLF_H1975_OSIvsControl    )))
geneDE_HCC827_OSIvsControl    <- as.data.frame(topTags(geneQLF_HCC827_OSIvsControl   , n = nrow(geneQLF_HCC827_OSIvsControl   )))
geneDE_HCC4006_OSIvsControl   <- as.data.frame(topTags(geneQLF_HCC4006_OSIvsControl  , n = nrow(geneQLF_HCC4006_OSIvsControl  )))



##Create a volcanoplot with enhanced_volcano_plot function filtering by logFC, pValue and RNA  
p1 <- enhanced_volcano_plot(geneDE_HCC4006_ERLOvsControl , "HCC4006"  ,  logFC_cutoff = 1,FDR_cutoff = 0.05, titlesize = 75,legendSize = 27)
p2 <- enhanced_volcano_plot(geneDE_HCC827_ERLOvsControl  , "HCC827"   ,  logFC_cutoff = 1,FDR_cutoff  = 0.05, titlesize = 75,legendSize = 27)
p3 <- enhanced_volcano_plot(geneDE_HCC827_OSIvsControl   , "HCC827"   ,  logFC_cutoff = 1,FDR_cutoff  = 0.05, titlesize = 75,legendSize = 27)
p5 <- enhanced_volcano_plot(geneDE_H1975_OSIvsControl    , "H1975"    ,  logFC_cutoff = 1,FDR_cutoff  = 0.05, titlesize = 75,legendSize = 27)
p4 <- enhanced_volcano_plot(geneDE_HCC4006_OSIvsControl  , "HCC4006"  ,  logFC_cutoff = 1,FDR_cutoff  = 0.05, titlesize = 75,legendSize = 27)

png("Volcano_plot_2.png", width =5300, height =1050)
grid.arrange(p1,p2, p3,p4,p5,ncol=5)
dev.off()


########Select Upregulated gnes in common in all tratments 
HCC827_ERLOvsControl_RNA_up   <- select_RNA(geneDE_HCC827_ERLOvsControl, FDR_value = 0.05,  logFC_Value = 1  )
HCC4006_ERLOvsControl_RNA_up  <- select_RNA(geneDE_HCC4006_ERLOvsControl,FDR_value = 0.05,  logFC_Value = 1  )
H1975_OSIvsControl_RNA_up     <- select_RNA(geneDE_H1975_OSIvsControl  , FDR_value = 0.05,  logFC_Value = 1  ) 
HCC827_OSIvsControl_RNA_up    <- select_RNA(geneDE_HCC827_OSIvsControl , FDR_value = 0.05,  logFC_Value = 1  ) 
HCC4006_OSIvsControl_RNA_up   <- select_RNA(geneDE_HCC4006_OSIvsControl, FDR_value = 0.05,  logFC_Value = 1  ) 



########Venn ERLO
venn_ERLO_upRNA <- list( HCC4006= HCC4006_ERLOvsControl_RNA_up$genes,
                         HCC827= HCC827_ERLOvsControl_RNA_up$genes)

#Create a list with upregulated RNA present in all groups
Venn_ERLO = Venn(venn_ERLO_upRNA)
upRNA_ERLO_list <- as.data.frame(overlap(Venn_ERLO, slice = "all"))
colnames(upRNA_ERLO_list) <- "Ensembl_ID_ERLO"

png ("Results/Venn_Diagram_ERLO.png", width = 1100, height = 850, units = "px")
ggvenn(Venn_ERLO,
       fill = c("red", "dodgerblue3")) + 
  theme_void() +
  theme(legend.position="none")
dev.off()

#obtain External gene name to RNA up-regulated
Biomart_common_RNA_ERLO <- get_gene_info(upRNA_ERLO_list )
write.table(file = "../Results/select_upRegulated_ERLO.csv", Biomart_common_RNA_ERLO, quote = F, sep= ",", row.names = F)

#######Venn OSI#####
venn_OSI_upRNA <-  list( HCC4006 = HCC4006_OSIvsControl_RNA_up$genes,  #list to obtain common RNA in all residual cells
                         HCC827  = HCC827_OSIvsControl_RNA_up$genes,
                         H1975   = H1975_OSIvsControl_RNA_up$genes)

#Create a list with upregulated RNA present in all groups
Venn_OSI = Venn(venn_OSI_upRNA)
upRNA_OSI_list <- as.data.frame(overlap(Venn_OSI, slice = "all"))
colnames(upRNA_OSI_list) <- "Ensembl_ID_OSI"

png ("/Results/Venn_Diagram_OSI.png", width = 1090, height = 1020, units = "px")
ggvenn(Venn_OSI,
       fill = c("red", "dodgerblue3", "deeppink")) + 
        theme_void() +
        theme(legend.position="none")
dev.off()


#obtain External gene name to RNA up-regulated
Biomart_common_RNA_OSI <- get_gene_info(upRNA_OSI_list )
write.table(file = "../Results/select_upRegulated_OSI.csv", Biomart_common_RNA_OSI, quote = F, sep= ",", row.names = F)

#######EGO Analysis####
ego_ERLO <- enrichGO(gene      = Biomart_common_RNA_ERLO$ensembl_gene_id,
                     OrgDb     = org.Hs.eg.db,
                     keyType   = 'ENSEMBL',
                     ont       = "ALL", pvalueCutoff = 0.5)

png("Results/Barplot_ERLO.png",width = 2300, height = 1850,res = 200)

barplot(ego_ERLO, font.size = 21.5,
                  title = "enrichGO erlotinib Overexpressed RNAs") +
                  theme(plot.title   = element_text(size = 23.5,face = "bold"),
                  legend.title = element_text(size = 23.5,face = "bold"), 
                  legend.text  =  element_text(size = 23.5,face = "bold"))
dev.off()


ego_Osi <- enrichGO(gene      = Biomart_common_RNA_OSI$ensembl_gene_id,
                    OrgDb     = org.Hs.eg.db,
                    keyType   = 'ENSEMBL',
                    ont       = "ALL", 
                    pvalueCutoff = 1)

png("/Results/Barplot_osi.png",width = 2300, height = 1850,res = 200)
 barplot(ego_Osi,
                       font.size = 21.5,
                       title = "enrichGO osimertinib Overexpressed RNAs") +
                 theme(plot.title   = element_text(size = 23.5,face = "bold"),
                 legend.title = element_text(size = 23.5,face = "bold"), 
                 legend.text  =  element_text(size = 23.5,face = "bold"))


dev.off()
#######Treeplot

c<- enrichplot::pairwise_termsim(ego_ERLO) 
dd<- enrichplot::pairwise_termsim(ego_Osi) 

treeplot(c)
treeplot(d)

sessionInfo()










########Select Upregulated gnes in common in all tratments 
HCC827_ERLOvsControl_RNA_down   <- select_RNA(geneDE_HCC827_ERLOvsControl, FDR_value = 0.05,  logFC_Value = -1  )
HCC4006_ERLOvsControl_RNA_down  <- select_RNA(geneDE_HCC4006_ERLOvsControl,FDR_value = 0.05,  logFC_Value = -1  )
H1975_OSIvsControl_RNA_down     <- select_RNA(geneDE_H1975_OSIvsControl  , FDR_value = 0.05,  logFC_Value = -1  ) 
HCC827_OSIvsControl_RNA_down    <- select_RNA(geneDE_HCC827_OSIvsControl , FDR_value = 0.05,  logFC_Value = -1  ) 
HCC4006_OSIvsControl_RNA_down   <- select_RNA(geneDE_HCC4006_OSIvsControl, FDR_value = 0.05,  logFC_Value = -1  ) 







########Venn ERLO
venn_ERLO_downRNA <- list( HCC4006= HCC4006_ERLOvsControl_RNA_down$genes,
                         HCC827= HCC827_ERLOvsControl_RNA_down$genes)

#Create a list with downregulated RNA present in all grodowns
Venn_ERLO_down = Venn(venn_ERLO_downRNA)
downRNA_ERLO_list <- as.data.frame(overlap(Venn_ERLO_down, slice = "all"))
colnames(downRNA_ERLO_list) <- "Ensembl_ID_ERLO"

#obtain External gene name to RNA up-regulated
Biomart_common_RNA_ERLO_down <- get_gene_info(downRNA_ERLO_list)
write.table(file = "Results/select_downRegulated_ERLO.csv", Biomart_common_RNA_ERLO_down, quote = F, sep= ",", row.names = F)


venn_OSI_downRNA <-  list( HCC4006 = HCC4006_OSIvsControl_RNA_down$genes,  #list to obtain common RNA in all residual cells
                           HCC827  = HCC827_OSIvsControl_RNA_down$genes,
                           H1975   = H1975_OSIvsControl_RNA_down$genes)

#Create a list with downregulated RNA present in all grodowns
Venn_OSI_down = Venn(venn_OSI_downRNA)
downRNA_OSI_list <- as.data.frame(overlap(Venn_OSI_down, slice = "all"))
colnames(downRNA_OSI_list) <- "Ensembl_ID_OSI"

#obtain External gene name to RNA down-regulated
Biomart_common_RNA_OSI_down <- get_gene_info(downRNA_OSI_list )
write.table(file = "Results/select_downRegulated_OSI.csv", Biomart_common_RNA_OSI_down, quote = F, sep= ",", row.names = F)

#######EGO Analysis####
ego_ERLO_down <- enrichGO(gene      = Biomart_common_RNA_ERLO_down$ensembl_gene_id,
                     OrgDb     = org.Hs.eg.db,
                     keyType   = 'ENSEMBL',
                     ont       = "all", pvalueCutoff = 0.05)
barplot(ego_ERLO_down,
  font.size = 21.5,
  title = "enrichGO erlotinib Overexpressed RNAs") +
  theme(plot.title   = element_text(size = 23.5,face = "bold"),
        legend.title = element_text(size = 23.5,face = "bold"), 
        legend.text  =  element_text(size = 23.5,face = "bold"))





ego_Osidown <- enrichGO(gene      = Biomart_common_RNA_OSI_dowmn$ensembl_gene_id,
                    OrgDb     = org.Hs.eg.db,
                    keyType   = 'ENSEMBL',
                    ont       = "ALL", 
                    pvalueCutoff = 0.05)

barplot(ego_Osidown,
        font.size = 21.5,
        title = "enrichGO osimertinib Overexpressed RNAs") +
  theme(plot.title   = element_text(size = 23.5,face = "bold"),
        legend.title = element_text(size = 23.5,face = "bold"), 
        legend.text  =  element_text(size = 13.5,face = "bold"))

png("Results/Barplot_down_osi.png",width = 2300, height = 1850,res = 200)
dd<- enrichplot::pairwise_termsim(ego_Osidown) 
treeplot(dd,   showCategory = 22) + ggtitle("enrichGO by downregulated genes in osimertinib residual cells") +
  theme(plot.title   = element_text(size = 23.5,face = "bold"), 
        legend.title = element_text(size = 13.5,face = "bold"))
dev.off()

sessionInfo()




