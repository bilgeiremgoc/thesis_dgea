gse <- getGEO("GSE120103")

gse1em <- getGEO("GSE120103", GSEMatrix = TRUE, AnnotGPL = TRUE)

gse_data <- gse1em[[1]]

exprs_data <- exprs(gse1em[[1]])

exprs_data <- log2(exprs_data + 1)

exprs(gse_data) <- exprs_data


boxplot(exprs_data, outline = FALSE, las=2, main="Before Normalization")

exprs_data <- normalizeBetweenArrays(exprs_data, method = "quantile")
boxplot(exprs_data, outline = FALSE, las=2, main="After Normalization")

metadata <- pData(gse_data)

group <- ifelse(grepl("endometriosis", metadata$title, ignore.case = TRUE), "control", "endometriosis")
group <- factor(group)
table(group)  # grup sayısını kontrol et


# Design matrisi oluştur
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)


fit <- lmFit(exprs_data, design)

# Grup kar????la??t??rmas??
contrast_matrix <- makeContrasts(endometriosis - control, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# Sonuçları alma

results <- topTable(fit2, adjust = "fdr", number = Inf)
head(results)

deg_results <- topTable(fit2, adjust.method = "BH", number = Inf)

# Filtrele: logFC > 1 ve adj.P.Val < 0.05
deg <- subset(deg_results, abs(logFC) > 1 & adj.P.Val < 0.05)
head(deg)

feature_data <- fData(gse1em[[1]])
colnames(feature_data)

deg$gene_symbol <- feature_data[rownames(deg), "Gene symbol"]
head(deg)

deg$gene_symbol[deg$gene_symbol == ""] <- NA

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
deg_rows <- rownames(deg_clean)[1:50]
exprs_subset <- exprs_data[deg_rows, ]
pheatmap(exprs_subset, scale = "row", main = "Heatmap of Top 50 DEGs")

write.csv(deg_clean, "deg_GSE120103.csv")


deg_filtered <- deg_clean %>%
  filter(adj.P.Val < 0.05 & abs(logFC) > 1)

up_genes <- deg_filtered %>% filter(logFC > 1)
down_genes <- deg_filtered %>% filter(logFC < -1)

cat("Upregulated:", nrow(up_genes), "\nDownregulated:", nrow(down_genes), "\n")










library(biomaRt)

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

gene_map <- getBM(
  filters = "external_gene_name",
  attributes = c("external_gene_name", "entrezgene_id"),
  values = deg_clean$gene_symbol,
  mart = mart
)

deg_clean <- merge(deg_clean, gene_map,
                   by.x = "gene_symbol", by.y = "external_gene_name",
                   all.x = TRUE)

deg_clean <- deg_clean[!is.na(deg_clean$entrezgene_id), ]

entrez_ids <- unique(deg_clean$entrezgene_id)

go_enrich <- enrichGO(
  gene          = entrez_ids,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

# Sonuçları görüntüle
head(go_enrich)

# Barplot (ilk 20 kategori)
barplot(go_enrich, showCategory = 20, title = "GO Biological Process Enrichment")

#kegg enrichment
kegg_enrich <- enrichKEGG(
  gene          = entrez_ids,
  organism      = "hsa",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05
)

# Entrez ID'den gene symbol'a çevir (okunabilirlik)
kegg_enrich <- setReadable(kegg_enrich, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

# İlk sonuçlara göz at
head(kegg_enrich)

# Barplot (ilk 20 yolak)
barplot(kegg_enrich, showCategory = 20, title = "KEGG Pathway Enrichment")

write_xlsx(as.data.frame(go_enrich), "GSE120103_GO_enrichment.xlsx")
write_xlsx(as.data.frame(kegg_enrich), "GSE120103_KEGG_enrichment.xlsx")

