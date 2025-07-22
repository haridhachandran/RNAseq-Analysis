# Load required packages ----
library(data.table)
library(DESeq2)
library(biomaRt)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ComplexHeatmap)
library(ggplot2)
library(pheatmap)

# ======= SET PARAMETERS =======

data_dir <- "../data/count_files"               # Folder containing sample count files
sample_metadata_file <- "../data/sample_metadata.csv"
results_dir <- "../results"
deg_dir <- file.path(results_dir, "DEGs")
figures_dir <- file.path(results_dir, "Figures")

dir.create(deg_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# =============================

# Step 1: Combine Count Files ----
file_list <- list.files(data_dir, pattern = "\\.txt$", full.names = TRUE)
data_frames <- lapply(file_list, fread)
combined_counts <- rbindlist(data_frames, idcol = "file", use.names = FALSE)
fwrite(combined_counts, file.path(results_dir, "combined_counts.tsv"), sep = "\t")

# Step 2: Load Sample Metadata ----
sample_table <- fread(sample_metadata_file)
rownames(sample_table) <- sample_table$sample_id

# Step 3: DESeq2 Object Creation ----
dds <- DESeqDataSetFromHTSeqCount(
  sampleTable = sample_table,
  directory = data_dir,
  design = ~ condition
)
dds <- DESeq(dds)

# Step 4: Extract Normalized Counts ----
normalized_counts <- as.data.frame(counts(dds, normalized = TRUE))
normalized_counts$ensembl_gene_id <- rownames(normalized_counts)

# Step 5: Gene Annotation ----
ensembl_ids <- rownames(normalized_counts)

# use biomaRt for annotation
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
annotation_info <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "description"),
  filters = "ensembl_gene_id",
  values = ensembl_ids,
  mart = ensembl
)

# Merge annotation
annotated_counts <- merge(normalized_counts, annotation_info, by = "ensembl_gene_id", all.x = TRUE)

# Step 6: DEG Analysis ----
res <- results(dds)
res_sig <- subset(res, padj < 0.05)

# DEG Subsets
up <- subset(res_sig, log2FoldChange > 0)
down <- subset(res_sig, log2FoldChange < 0)
res_strict <- subset(res_sig, abs(log2FoldChange) > 1)

# Save DEGs
write.csv(res_strict, file = file.path(deg_dir, "significant_DEGs_lfc1.csv"))
write.csv(up, file = file.path(deg_dir, "upregulated_genes.csv"))
write.csv(down, file = file.path(deg_dir, "downregulated_genes.csv"))

# Step 7: Plots ----

# 1. MA Plot
pdf(file.path(figures_dir, "MA_plot.pdf"))
plotMA(res, ylim = c(-5, 5))
dev.off()

# 2. Volcano Plot
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
res_df$Significant <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, "Yes", "No")

volcano <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = Significant)) +
  scale_color_manual(values = c("gray", "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot", x = "Log2FC", y = "-log10(padj)")

ggsave(file.path(figures_dir, "Volcano_plot.png"), volcano, width = 7, height = 5)

# 3. PCA Plot
rld <- rlog(dds, blind = TRUE)
pcaData <- plotPCA(rld, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
p <- ggplot(pcaData, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA of Samples") +
  theme_minimal()

ggsave(file.path(figures_dir, "PCA_plot.png"), p, width = 7, height = 5)

# 4. Heatmap of Top Variable Genes
select_genes <- head(order(rowVars(assay(rld)), decreasing = TRUE), 50)
mat <- assay(rld)[select_genes, ]
mat <- t(scale(t(mat)))
annot <- as.data.frame(colData(rld)[, "condition", drop = FALSE])
pheatmap(mat, annotation_col = annot, filename = file.path(figures_dir, "Heatmap_top50.png"))

