# Load the DESeq2 library
library(DESeq2)
library(tidyverse)
library(EnhancedVolcano)
library(gprofiler2)

# Load the counts matrix and metadata table from a CSV file
counts_matrix <- read.csv('counts_chinmo.csv', row.names=1)
metadata <- read.csv("chinmo_metadata.csv", row.names=1)

# Create a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = counts_matrix,
                              colData = metadata,
                              design = ~condition,
                              tidy = TRUE)

# Perform differential expression analysis
dds <- DESeq(dds)

# Get the results table
results <- results(dds, contrast = c("condition", "experimental", "control"), tidy = TRUE)

# Create a volcano plot to visualize differential expression
png(filename = "chinmo_RNAi_volcano_plot.png", height = 1080, width = 1100)
plot <- EnhancedVolcano(results,
                        lab = results$row,
                        x = 'log2FoldChange',
                        y = 'pvalue',
                        title = "Chinmo RNAi",
                        pointSize = 1.5,
                        labSize = 6.0,
                        pCutoff = 0.05,
                        FCcutoff = 0.5,
                        colAlpha = 1,
                        legendLabels = c('NotSig', 'Log2FC', 'pvalue',
                                         'pvalue & Log2FC'),
                        legendPosition = 'right',
                        legendLabSize = 16,
                        legendIconSize = 6.0) +
  ggplot2::coord_cartesian(xlim = c(-10, 10))
print(plot)
dev.off()

# Get lists of significantly upregulated and downregulated genes
significant_results <- results[results$padj<0.05,]

upreg_data <- significant_results %>%
  filter(log2FoldChange > 0) %>%
  arrange(pvalue)
upregulated_genes <- upreg_data$row

downreg_data <- significant_results %>%
  filter(log2FoldChange < 0) %>%
  arrange(pvalue)
downregulated_genes <- downreg_data$row

# perform functional profiling
gostres_up <- gost(upregulated_genes, ordered_query=TRUE, organism="dmelanogaster", correction_method="gSCS", 
                   domain_scope="annotated", sources="GO:BP") 
gost_results_upreg <- gostres_up$result[,c("term_name","p_value")]
write.csv(gost_results_upreg, "chinmo_gost_results_upreg.csv")

gostres_down <- gost(downregulated_genes, ordered_query=TRUE, organism="dmelanogaster", correction_method="gSCS", 
                     domain_scope="annotated", sources="GO:BP")
gost_results_downreg <- gostres_down$result[,c("term_name", "p_value")]
write.csv(gost_results_downreg, "chinmo_gost_results_downreg.csv")


