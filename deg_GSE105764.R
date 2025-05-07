library(GEOquery)
library(DESeq2)
library(ggplot2)
library(pheatmap)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GEOquery")

library(GEOquery)

url <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE105nnn/GSE105764/suppl/GSE105764_count.txt.gz"
download.file(url, destfile = "GSE105764_count.txt.gz")


counts <- read.table("GSE105764_count.txt.gz", header = TRUE, row.names = 1)


print(dim(counts))  # Gen sayısı x örnek sayısı
print(head(counts))  # İlk birkaç satır

pheno <- pData(data)
group <- ifelse(grepl("case", pheno$characteristics_ch1, ignore.case = TRUE), "control", "case")
group <- factor(group)
table(group)  


dge <- DGEList(counts = counts, group = group)

# Düşük ifadeli genleri filtrele
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes=FALSE]

# Veriyi normalize et (TMM yöntemi)
dge <- calcNormFactors(dge)

# Deney tasarımı
design <- model.matrix(~group)
dge <- estimateDisp(dge, design)

# Farklı ifadelenen genleri bul
fit <- glmFit(dge, design)
lrt <- glmLRT(fit, coef=2)

# Sonuçları al ve anlamlı genleri filtrele
deg_results <- topTags(lrt, n = Inf)
deg_results <- deg_results$table

gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = rownames(deg_results),
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

deg_results$GeneSymbol <- gene_symbols

head(deg_results)

# Anlamlı genleri filtrele (FDR < 0.05 ve mutlak logFC > 1)
deg_filtered <- deg_results[deg_results$FDR < 0.05 & abs(deg_results$logFC) > 1, ]

# Sıralı listeyi yazdır
head(deg_filtered)

# Anlamlı gen sayısı
cat("Anlamlı gen sayısı:", nrow(deg_filtered), "\n")


write_xlsx(deg_filtered, "deg_GSE105764.xlsx")

deg_filtered <- deg_results %>%
  filter(FDR < 0.05 & abs(logFC) > 1)

# Kaç tane gen kaldı?
nrow(deg_filtered)

# İlk birkaç tanesini gör
head(deg_filtered)

up_genes <- deg_filtered %>% filter(logFC > 1)
down_genes <- deg_filtered %>% filter(logFC < -1)



# Volcano plot
deg_results$significant <- ifelse(deg_results$FDR < 0.05 & abs(deg_results$logFC) > 1, "Significant", "Not Significant")

ggplot(deg_results, aes(x = logFC, y = -log10(FDR), color = significant)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("grey", "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot - GSE105764",
       x = "Log2 Fold Change",
       y = "-Log10 FDR") +
  theme(legend.position = "top")



library(EnhancedVolcano)

EnhancedVolcano(deg_results,
                lab = rownames(deg_results),
                x = 'logFC',
                y = 'PValue',
                pCutoff = 0.05,
                FCcutoff = 1,
                title = "Volcano Plot for GSE105764",
                subtitle = "Differential Expression Analysis")



library(pheatmap)

# En anlamlı genlerden 20'sini seçme
top_genes <- rownames(deg_results)[1:20]
heatmap_data <- counts[top_genes, ]
heatmap_data <- log2(heatmap_data + 1)  # Log dönüşümü

pheatmap(heatmap_data, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         scale = "row", 
         show_rownames = TRUE,
         main = "Top 20 DEGs")


gene_list <- rownames(deg_results)
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
barplot(go_enrich, showCategory = 20, title = "GO BP Enrichment")

write_xlsx(as.data.frame(go_enrich), "GSE105764_GO_enrichment.xlsx")

# KEGG analizi
kegg_enrich <- enrichKEGG(gene = entrez_ids,
                          organism = "hsa",
                          pAdjustMethod = "BH",
                          qvalueCutoff = 0.05)


# Barplot
dotplot(kegg_enrich, showCategory=10, title="Top 10 KEGG Enriched Pathways")

write_xlsx(as.data.frame(kegg_enrich), "GSE105764_KEGG_enrichment.xlsx")
