# Load the DESeq2 library
library(DESeq2)

# Load the counts matrix from a CSV file
counts_matrix <- read.csv('counts_matrix.csv', row.names=1)

# Create a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = counts_matrix,
                              colData = NULL,
                              design = ~ condition)

# Perform differential expression analysis
dds <- DESeq(dds)

# Get the results table
results <- results(dds)

# Create a volcano plot to visualize differential expression
volcano_plot <- function(results) {
  plot(results$log2FoldChange, -log10(results$pvalue),
       xlab = "Log2 Fold Change", ylab = "-log10 p-value",
       main = "Volcano Plot of Differential Expression",
       xlim = c(-5, 5))
  # Highlight significant genes
  sig_genes <- subset(results, padj < 0.05)
  points(sig_genes$log2FoldChange, -log10(sig_genes$pvalue), col = "red", pch = 16)
  legend("topright", legend = c("Significant Genes"), col = c("red"), pch = c(16), cex = 0.8)
}

# Plot the volcano plot
volcano_plot(results)