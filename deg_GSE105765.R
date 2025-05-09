
library(DESeq2)
library(ggplot2)
library(pheatmap)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GEOquery")

library(GEOquery)

gse <- getGEO("GSE105765")
getGEOSuppFiles("GSE105765")
gse1em <- getGEO("GSE105765", GSEMatrix = TRUE)

files <- list.files("GSE105765", pattern = "txt.gz", full.names = TRUE)

# Her bir örneği oku
counts_list <- lapply(files, function(f) {
  dat <- read.table(gzfile(f), header = TRUE, sep = "\t", row.names = 1)
  dat[, 1, drop = FALSE]  # ilk sütun count gibi görünüyorsa
})

# Ortak satır isimlerine göre birleştir
counts <- do.call(cbind, counts_list)
colnames(counts) <- gsub("_.*", "", basename(files))  # Kolon isimleri


exprs_data <- exprs(gse[[1]])
meta <- pData(gse[[1]])


counts <- read.table("C:/Users/bilge/Documents/thesis_dgea/GSE105765_count.txt.gz", 
                     header = TRUE, 
                     sep = "\t", 
                     row.names = 1, 
                     check.names = FALSE)
dim(counts)  # Satır ve sütun sayısını kontrol et
colnames(counts)

group <- ifelse(grepl("^EC", colnames(counts)), "control", "case")
group <- factor(group)
table(group)  


# colData oluşturma
colData <- data.frame(row.names = colnames(counts),
                      condition = factor(group, levels = c("control", "case")))

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = colData,
                              design = ~ condition)



dds <- DESeq(dds)

res <- results(dds, contrast = c("condition", "control", "case"))

deg <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
head(deg)

result <- miRNA_NameToAccession(rownames(deg), version = "v22")
head(result)

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

data <- read.csv("deg_GSE105765.csv")

data_cleaned <- na.omit(data)
head(data_cleaned)

write_xlsx(data_cleaned, "deg_GSE105765.xlsx")



rownames(data_cleaned$symbol)
gene_list <- data_cleaned$X
head(gene_list)
gene_df <- bitr(gene_list, fromType="ENSEMBL", toType=c("SYMBOL", "ENTREZID"), OrgDb="org.Hs.eg.db")
head(gene_df)

# GO analizi
go_results <- enrichGO(gene = gene_df$ENTREZID,
                       OrgDb = org.Hs.eg.db,
                       ont = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05,
                       readable = TRUE)
head(go_results)

#Barplot
barplot(go_results, showCategory = 20, title = "GO BP Enrichment")

write_xlsx(as.data.frame(go_results), "GSE105765_GO_enrichment.xlsx")

# KEGG analizi
kegg_enrich <- enrichKEGG(gene = gene_df$ENTREZID,
                          organism = "hsa",
                          pAdjustMethod = "BH",
                          qvalueCutoff = 0.05)


# Barplot
dotplot(kegg_enrich, showCategory=10, title="Top 10 KEGG Enriched Pathways")

write_xlsx(as.data.frame(kegg_enrich), "GSE105765_KEGG_enrichment.xlsx")


deg_filtered <- data_cleaned %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1)

# Kaç tane gen kaldı?
nrow(deg_filtered)

# İlk birkaç tanesini gör
head(deg_filtered)

up_genes <- deg_filtered %>% filter(log2FoldChange > 1)
down_genes <- deg_filtered %>% filter(log2FoldChange < -1)
