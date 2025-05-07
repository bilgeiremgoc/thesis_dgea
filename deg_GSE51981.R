library(GEOquery)
library(limma)
library(annotate)
library(ggplot2)
library(pheatmap)
library(hgu133a.db)
library(dplyr)


gse <- getGEO("GSE51981")

gse1em <- getGEO("GSE51981", GSEMatrix = TRUE, AnnotGPL = TRUE)

exprs_data <- exprs(gse1em[[1]])

if (max(exprs_data, na.rm = TRUE) > 100) {
  exprs_data <- log2(exprs_data + 1)
}


boxplot(exprs_data, outline = FALSE, las=2, main="Before Normalization")

exprs_data <- normalizeBetweenArrays(exprs_data)
boxplot(exprs_data, outline = FALSE, las=2, main="After Normalization")

gse_data<- gse[[1]]

metadata <- pData(gse_data)
group <- ifelse(metadata$`endometriosis/no endometriosis:ch1` == "Endometriosis", "case", "control")
group <- factor(group)
table(group)


# Design matrisi oluştur
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)


#Normalizasyon var mı kontrol etme ve normalize etme
boxplot(exprs_data, outline = FALSE, las=2, main="Before Normalization")

exprs_data <- normalizeBetweenArrays(exprs_data)
boxplot(exprs_data, outline = FALSE, las=2, main="After Normalization")

fit <- lmFit(exprs_data, design)
colnames(design)
# Grup kar????la??t??rmas??
contrast_matrix <- makeContrasts(case - control, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# Sonu??lar?? al

results <- topTable(fit2, adjust = "fdr", number = Inf)

deg <- subset(results, abs(logFC) > 1 & adj.P.Val < 0.05)
head(deg)

# Volkan Plotu
results$threshold <- as.factor(abs(results$logFC) > 1 & results$adj.P.Val < 0.05)

ggplot(results, aes(x = logFC, y = -log10(adj.P.Val), color = threshold)) +
  geom_point(alpha = 0.8) +
  theme_minimal() +
  scale_color_manual(values = c("grey", "red")) +
  labs(title = "Volcano Plot")

deg_rows <- rownames(deg)[1:50]

# Heatmap (ilk 50 DEG)
pheatmap(exprs_data[deg_rows, ], scale = "row", main = "Heatmap of Top 50 DEGs")

deg$GeneSymbol <- getSYMBOL(rownames(deg), "hgu133a.db")
head(deg)

deg_clean <- deg[!is.na(deg$GeneSymbol), ]
deg_clean


write.csv(deg_clean, "deg_GSE51981.csv")


deg_filtered <- deg_clean %>%
  filter(adj.P.Val < 0.05 & abs(logFC) > 1)

up_genes <- deg_filtered %>% filter(logFC > 1)
down_genes <- deg_filtered %>% filter(logFC < -1)

cat("Upregulated:", nrow(up_genes), "\nDownregulated:", nrow(down_genes), "\n")










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
entrez_loose <- na.omit(as.character(deg_loose$ENTREZID))

# enrichKEGG fonksiyonunu çalıştır (burası eksik şu anda)
kegg_enrich <- enrichKEGG(gene = entrez_loose,
                          organism = "hsa",
                          pvalueCutoff = 0.1)

# Gen adlarına çevir (okunabilir hale getir)
kegg_enrich <- setReadable(kegg_enrich, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

# Artık çıktı bakılabilir
head(kegg_enrich)

# Barplot
barplot(kegg_enrich, showCategory = 20, title = "KEGG Pathway Enrichment")

write_xlsx(as.data.frame(go_enrich), "GSE51981_GO_enrichment.xlsx")
write_xlsx(as.data.frame(kegg_enrich), "GSE51981_KEGG_enrichment.xlsx")

