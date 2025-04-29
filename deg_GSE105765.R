
library(DESeq2)
library(ggplot2)
library(pheatmap)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GEOquery")

library(GEOquery)

gse <- getGEO("GSE105765")

gse1em <- getGEO("GSE105765", GSEMatrix = TRUE, AnnotGPL = TRUE)

meta <- pData(gse[[1]])

# Örnek grup bilgisi
group <- meta$`disease state`  # veya uygun olan sütun
colData <- data.frame(row.names = colnames(counts),
                      condition = factor(group))

# 1. Sütunu faktör olarak ayarla ve referans seviyesini 'control' yap
colData$condition <- factor(colData$condition, levels = c("control", "case"))

# 2. Ardından DESeqDataSet oluştur
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = colData,
                              design = ~ condition)



dds <- DESeq(dds)

res <- results(dds, contrast = c("condition", "control", "case"))

deg <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
head(deg)

library(org.Hs.eg.db)
deg$symbol <- mapIds(org.Hs.eg.db,
                     keys = rownames(deg),
                     column = "SYMBOL",
                     keytype = "ENSEMBL",
                     multiVals = "first")


# 6. Volcano plot
res$threshold <- as.factor(res$padj < 0.05 & abs(res$log2FoldChange) > 1)

volcano_df <- na.omit(res)  # NA olan satırları çıkar

volcano_df$threshold <- as.factor(volcano_df$padj < 0.05 & abs(volcano_df$log2FoldChange) > 1)

ggplot(volcano_df, aes(x = log2FoldChange, y = -log10(padj), color = threshold)) +
  geom_point(alpha = 0.8) +
  theme_minimal() +
  scale_color_manual(values = c("grey", "red")) +
  labs(title = "Volcano Plot")




# 7. Heatmap
vsd <- vst(dds)
pheatmap(assay(vsd)[rownames(deg)[1:50], ], scale = "row",
         annotation_col = colData)

write.csv(deg, "deg_GSE105765.csv" )

deg
