# scripts/deseq2_analysis.R

# Load packages
library(DESeq2)
library(EnhancedVolcano)
library(pheatmap)
library(ggplot2)

# Load the count data
countData <- read.table("data/gene_counts.txt", header = TRUE, row.names = 1, sep = "\t")

# Use columns 6 to 11 which are 3 control and 3 LPS samples
countData <- countData[, 6:11]

# Create metadata for the samples
condition <- factor(c("control", "control", "control", "LPS", "LPS", "LPS"))
colData <- data.frame(row.names = colnames(countData), condition)

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = round(countData),
                              colData = colData,
                              design = ~ condition)

# Filter low counts
dds <- dds[rowSums(counts(dds)) > 10, ]

# Run DESeq2
dds <- DESeq(dds)
res <- results(dds)

# Save results
write.csv(as.data.frame(res), file = "results/deseq2_results.csv")

# Volcano Plot
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'LPS vs Control')

# PCA Plot
vsd <- vst(dds)
pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal()

# Heatmap of top 30 genes
select <- order(rowMeans(counts(dds, normalized = TRUE)), decreasing = TRUE)[1:30]
pheatmap(assay(vsd)[select,], cluster_rows = TRUE, cluster_cols = TRUE,
         annotation_col = colData)

# Save Volcano Plot - High-resolution PNG (300 DPI, 8 x 8 inches)
png("results/volcano_plot.png", width = 8, height = 8, units = "in", res = 300)
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'LPS vs Control')
dev.off()

# PCA Plot
png("results/pca_plot.png", width = 8, height = 6, units = "in", res = 300)
ggplot(pcaData, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal()
dev.off()

# Heatmap visualization
# Set annotation color for control and LPS
ann_colors <- list(
  condition = c(control = "#00BFC4", LPS = "#F8766D")
)

# Save high-resolution heatmap
png("results/heatmap.png", width = 10, height = 8, units = "in", res = 300)

pheatmap(
  assay(vsd)[select, ],
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = colData,
  show_rownames = TRUE,           # Show gene IDs
  fontsize_row = 7,               # Reduce row font size
  fontsize_col = 12,              # Bigger sample labels
  angle_col = 45,                 # Tilt sample labels
  main = "Top 30 Expressed Genes",
  annotation_colors = ann_colors
)

dev.off()


