gse <- getGEO("GSE23339")

gse1em <- getGEO("GSE23339", GSEMatrix = TRUE, AnnotGPL = TRUE)

exprs_data <- exprs(gse1em[[1]])

exprs_data <- log2(exprs_data + 1)


boxplot(exprs_data, outline = FALSE, las=2, main="Before Normalization")

#Control ve Hastalık gruplarını oluşturma
gse_data<- gse[[1]]

exprs(gse_data) <- exprs_data

metadata <- pData(gse_data)
group <- ifelse(grepl("control", metadata$title, ignore.case = TRUE), "control", "endometriosis")
group <- factor(group)
table(group)  # grup sayısını kontrol et

design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

fit <- lmFit(exprs_data, design)
contrast_matrix <- makeContrasts(endometriosis - control, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

results <- topTable(fit2, adjust = "fdr", number = Inf)

deg <- subset(results, abs(logFC) > 1 & P.Value < 0.05)
head(deg)

# Volkan Plotu
results$threshold <- as.factor(abs(results$logFC) > 1 & results$adj.P.Val < 0.05)

ggplot(results, aes(x = logFC, y = -log10(adj.P.Val), color = threshold)) +
  geom_point(alpha = 0.8) +
  theme_minimal() +
  scale_color_manual(values = c("grey", "red")) +
  labs(title = "Volcano Plot")

# Heatmap (ilk 50 DEG)
pheatmap(deg_filtered[rownames(deg)[1:50], ], scale = "row")

feature_data <- fData(gse1em[[1]])

head(feature_data[ , 1:5])
results$GeneSymbol <- feature_data[rownames(results), "Gene symbol"]

deg <- results[!is.na(results$GeneSymbol) & results$GeneSymbol != "", ]
head(deg)


write_xlsx(deg, "deg_GSE23339.xlsx")

deg_filtered <- deg %>%
  filter(adj.P.Val < 0.05 & abs(logFC) > 1)

# Kaç tane gen kaldı?
nrow(deg_filtered)

# İlk birkaç tanesini gör
head(deg_filtered)

up_genes <- deg_filtered %>% filter(logFC > 1)
down_genes <- deg_filtered %>% filter(logFC < -1)


if (!requireNamespace("illuminaHumanv4.db", quietly = TRUE)) {
  BiocManager::install("illuminaHumanv4.db")
}
library(illuminaHumanv4.db)








deg_filtered$ENTREZID <- mapIds(illuminaHumanv4.db,
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

write_xlsx(as.data.frame(go_enrich), "deg_GSE23339_GO_enrichment.xlsx")
write_xlsx(as.data.frame(kegg_enrich), "deg_GSE23339_KEGG_enrichment.xlsx")


