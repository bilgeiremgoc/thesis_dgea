# Gerekli kütüphaneler
library(ggplot2)
library(pheatmap)
library(dplyr)

gse <- getGEO("GSE7307")

gse1em <- getGEO("GSE7307", GSEMatrix = TRUE, AnnotGPL = TRUE)

exprs_data <- exprs(gse1em[[1]])


if (max(exprs_data, na.rm = TRUE) > 100) {
  exprs_data <- log2(exprs_data + 1)
}


boxplot(exprs_data, outline = FALSE, las=2, main="Before Normalization")

exprs_data <- normalizeBetweenArrays(exprs_data)
boxplot(exprs_data, outline = FALSE, las=2, main="After Normalization")

gse_data<- gse[[1]]

metadata <- pData(gse_data)


head(metadata$characteristics_ch1)

# Endometrium doku olanları seç
endometrium_samples <- grepl("endometrium", metadata$characteristics_ch1, ignore.case = TRUE)

# Veriyi filtrele
exprs_data <- exprs_data[, endometrium_samples]
metadata <- metadata[endometrium_samples, ]

group <- ifelse(grepl("normal", metadata$`Disease type:ch1`, ignore.case = TRUE), "endometriosis", "normal")
group <- factor(group)
table(group)

design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

fit <- lmFit(exprs_data, design)
contrast_matrix <- makeContrasts(endometriosis - normal, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)


results <- topTable(fit2, adjust = "fdr", number = Inf)
head(results)

deg_results <- topTable(fit2, adjust.method = "BH", number = Inf)
# Filtrele: logFC > 1 ve adj.P.Val < 0.05
deg <- subset(deg_results, abs(logFC) > 1 & adj.P.Val < 0.05)
head(deg)

deg$gene_symbol <- getSYMBOL(rownames(deg), "hgu133a.db")
head(deg)

deg_clean <- deg[!is.na(deg$gene_symbol), ]
deg_clean

# Volkan Plotu
results$threshold <- as.factor(abs(results$logFC) > 1 & results$adj.P.Val < 0.05)

ggplot(results, aes(x = logFC, y = -log10(adj.P.Val), color = threshold)) +
  geom_point(alpha = 0.8) +
  theme_minimal() +
  scale_color_manual(values = c("grey", "red")) +
  labs(title = "Volcano Plot")


# Heatmap (ilk 50 DEG)
pheatmap(exprs_data[rownames(deg_clean)[1:50], ], scale = "row")
#Annotation sorunu yaşandığı için çıkmadı ama gen listesi var



write_xlsx(deg_clean, "deg_GSE7307.xlsx")

deg_filtered <- deg_clean %>%
  filter(adj.P.Val < 0.05 & abs(logFC) > 1)

# Kaç tane gen kaldı?
nrow(deg_filtered)

# İlk birkaç tanesini gör
head(deg_filtered)

up_genes <- deg_filtered %>% filter(logFC > 1)
down_genes <- deg_filtered %>% filter(logFC < -1)







deg_filtered$ENTREZID <- mapIds(hgu133a.db,
                                keys = rownames(deg_filtered),
                                column = "ENTREZID",
                                keytype = "PROBEID",
                                multiVals = "first")



# NA olmayanları al
entrez_ids <- na.omit(deg_filtered$ENTREZID)

go_enrich <- enrichGO(gene = entrez_ids,
                      OrgDb = org.Hs.eg.db,
                      keyType = "ENTREZID",
                      ont = "BP",
                      pAdjustMethod = "BH",
                      qvalueCutoff = 0.05,
                      readable = TRUE)

head(go_enrich)

#Barplot
barplot(go_enrich, showCategory = 20, title = "GO BP Enrichment")

#kegg enrichment

kegg_enrich <- enrichKEGG(gene = entrez_ids,
                          organism = "hsa",
                          pAdjustMethod = "BH",
                          qvalueCutoff = 0.05)

# Entrez ID'den gen adlarına dönüştür (okunabilirlik)
kegg_enrich <- setReadable(kegg_enrich, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

# İlk sonuçlara bak
head(kegg_enrich)

# Barplot
barplot(kegg_enrich, showCategory = 20, title = "KEGG Pathway Enrichment")

write_xlsx(as.data.frame(go_enrich), "GSE7307_GO_enrichment.xlsx")
write_xlsx(as.data.frame(kegg_enrich), "GSE7307_KEGG_enrichment.xlsx")
