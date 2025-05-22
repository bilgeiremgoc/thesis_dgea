gse <- getGEO("GSE87809", GSEMatrix = TRUE)

gse <- gse[[1]]


countData <- read.delim("GSE87809_raw_counts_GRCh38.p13_NCBI.tsv.gz", row.names = 1)
head(countData)

metadata <- pData(gse)
head(metadata)

group_list <- ifelse(grepl("control", metadata$characteristics_ch1.1), "control", "endometriosis")
group_list <- factor(group_list)
table(group_list)

# Tasarım matrisi
design <- model.matrix(~ 0 + group_list)
colnames(design) <- levels(group_list)

dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = data.frame(condition = group_list),
                              design = ~ condition)

dds <- DESeq(dds)
res <- results(dds)
resOrdered <- res[order(res$padj), ]

resSig <- subset(resOrdered, padj < 0.05 & abs(log2FoldChange) > 1)

# DEG listesi dosyaya kaydedelim
write.csv(as.data.frame(resSig), "GSE87809_DEG_list.csv")

resOrdered$threshold <- ifelse(resOrdered$padj < 0.05 & abs(resOrdered$log2FoldChange) > 1, "Significant", "Not Significant")

gene_names <- featureNames(gse)  # GEO veri setindeki gen isimlerini al
names(gene_names) <- gene_names

resSig$gene_id <- rownames(resSig)
resSig$gene_name <- gene_names[resSig$gene_id]

# Volcano plot oluştur
ggplot(resOrdered, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = threshold)) +
  theme_minimal() +
  labs(title = "Volcano Plot", x = "log2 Fold Change", y = "-log10 p-value") +
  scale_color_manual(values = c("blue", "red"))


sig_genes <- rownames(resSig)
exprSig <- countData[sig_genes, ]

exprSig_log2 <- log2(exprSig + 1)

pheatmap(exprSig_log2,
         scale = "row",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         clustering_method = "complete",
         show_rownames = FALSE,
         main = "Heatmap of Significant DEGs")



# resSig'i data frame'e dönüştürüp, Excel dosyasına kaydedelim
resSig_df <- as.data.frame(resSig)

# Excel dosyasına yazdıralım
write_xlsx(resSig_df, "GSE87809_DEG_with_genenames.xlsx")


resSig$clean_id <- gsub("\\..*", "", rownames(resSig))

library(org.Hs.eg.db)

resSig$gene_name <- mapIds(org.Hs.eg.db,
                           keys = resSig$clean_id,
                           column = "SYMBOL",
                           keytype = "ENTREZID",
                           multiVals = "first")

write_xlsx(as.data.frame(resSig), "GSE87809_DEG_with_gene_names.xlsx")



# Gen tanımını ekle
resSig$gene_description <- mapIds(org.Hs.eg.db,
                                  keys = resSig$clean_id,
                                  column = "GENENAME",
                                  keytype = "ENTREZID",
                                  multiVals = "first")

# Göz at
head(resSig[, c("clean_id", "gene_name", "gene_description")])

write_xlsx(as.data.frame(resSig), "GSE87809_DEG_with_gene_names_and_description.xlsx")



#GENE ONTOLOGY AND KEGG DATABASE RESULTS

if (!require("clusterProfiler")) BiocManager::install("clusterProfiler")
library(clusterProfiler)
library(org.Hs.eg.db)

# Entrez ID'leri al
entrez_ids <- resSig$clean_id[!is.na(resSig$clean_id)]

# GO analizi: Biyolojik süreç (BP)
go_enrich <- enrichGO(gene = entrez_ids,
                      OrgDb = org.Hs.eg.db,
                      keyType = "ENTREZID",
                      ont = "BP",
                      pAdjustMethod = "BH",
                      qvalueCutoff = 0.05,
                      readable = TRUE)

head(go_enrich)

barplot(go_enrich, showCategory = 20, title = "Top 20 GO Biological Processes")

# KEGG analizi
kegg_enrich <- enrichKEGG(gene = entrez_ids,
                          organism = 'hsa',
                          pAdjustMethod = "BH",
                          qvalueCutoff = 0.05)

# Sonuçları daha okunabilir hâle getir
kegg_enrich <- setReadable(kegg_enrich, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

# İlk birkaç sonucu gör
head(kegg_enrich)

barplot(kegg_enrich, showCategory = 20, title = "Top 20 KEGG Pathways")

write_xlsx(as.data.frame(go_enrich), "GSE87809_GO_enrichment_results.xlsx")
write_xlsx(as.data.frame(kegg_enrich), "GSE87809_KEGG_enrichment_results.xlsx")



