library(GEOquery)
library(DESeq2)
library(ggplot2)
library(pheatmap)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GEOquery")

library(GEOquery)

gse <- getGEO("GSE105764")


#database error verdi tekrar dene

getGEOSuppFiles("GSE105764")

# Sayım matrisi oku
counts <- read.table("C:/Users/bilge/Documents/thesis_dgea/GSE105764/GSE105764_count.txt.gz",
                     header = TRUE, row.names = 1, sep = "\t")
head(counts)


group <- factor(c(rep("endometriosis", 8), rep("control", 8)))


colData <- data.frame(
  row.names = colnames(counts),
  condition = group
)

counts <- read.table("C:/Users/bilge/Documents/thesis_dgea/GSE105764/GSE105764_count.txt.gz",
                     header = TRUE, row.names = 1, sep = "\t")

# 2. Annotation kolonlarını kaldır (eğer varsa)
# Eğer sadece sayım varsa bu adım atlanabilir:
# counts <- counts[, 7:ncol(counts)]  # Gerekirse

# 3. Grup bilgisi
group <- factor(c(rep("endometriosis", 8), rep("control", 8)))
colData <- data.frame(condition = group)
rownames(colData) <- colnames(counts)

# 4. DESeq2 analiz
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = colData,
                              design = ~ condition)
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "endometriosis", "control"))

# 5. DEGs
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


write.csv(deg, "deg_GSE105764.csv" )

deg
