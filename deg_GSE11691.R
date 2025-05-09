install.packages("BiocManager")
BiocManager::install("GEOquery")
library(GEOquery)

BiocManager::install(c("GEOquery", "limma", "annotate", "hgu133a.db","clusterProfiler"))
install.packages("ggplot2")
install.packages("pheatmap")


library(GEOquery)
library(limma)
library(annotate)
library(hgu133a.db)  # Platform GPL96 için uygun
library(ggplot2)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(DOSE)

#Dataseti elde etme ve matrix oluşturma
gse <- getGEO("GSE11691")
gse1em <- getGEO("GSE11691", GSEMatrix = TRUE, AnnotGPL = TRUE)

gse_data <- gse1em[[1]]

exprs_data <- exprs(gse1em[[1]])

exprs_data <- log2(exprs_data + 1)

exprs(gse_data) <- exprs_data


metadata <- pData(gse_data)
group <- ifelse(grepl("endometriosis", metadata$title, ignore.case = TRUE), "control", "endometriosis")
group <- factor(group)
table(group)  # grup sayısını kontrol et



# Design matrisi oluştur
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)


#Normalizasyon var mı kontrol etme ve normalize etme
boxplot(exprs_data, outline = FALSE, las=2, main="Before Normalization")

exprs_data <- normalizeBetweenArrays(exprs_data)
boxplot(exprs_data, outline = FALSE, las=2, main="After Normalization")

fit <- lmFit(exprs_data, design)

# Grup karşılaştırması
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

deg$gene_symbol <- getSYMBOL(rownames(deg), "hgu133a.db")
head(deg)


# Volkan Plotu

deg_results$threshold <- as.factor(abs(deg_results$logFC) > 1 & deg_results$adj.P.Val < 0.05)

ggplot(deg_results, aes(x = logFC, y = -log10(adj.P.Val), color = threshold)) +
  geom_point(alpha = 0.8, size = 1.5) +
  scale_color_manual(values = c("grey", "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot",
       x = "logFC",
       y = "-log10(adjusted P-value)",
       color = "Significant")

# Heatmap (ilk 50 DEG)
pheatmap(exprs_data[rownames(deg)[1:50], ], scale = "row")


write_xlsx(deg, "deg_GSE11691.xlsx")

deg_filtered <- deg_results %>%
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


kegg_enrich <- setReadable(kegg_enrich, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

head(kegg_enrich)

# Barplot
barplot(kegg_enrich, showCategory = 20, title = "KEGG Pathway Enrichment")

write_xlsx(as.data.frame(go_enrich), "GSE11691_GO_enrichment.xlsx")
write_xlsx(as.data.frame(kegg_enrich), "GSE11691_KEGG_enrichment.xlsx")



