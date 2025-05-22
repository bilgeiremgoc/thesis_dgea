gse <- getGEO("GSE120103")

gse1em <- getGEO("GSE120103", GSEMatrix = TRUE, AnnotGPL = TRUE)

gse_data <- gse1em[[1]]

exprs_data <- exprs(gse_data)

boxplot(exprs_data, outline = FALSE, las=2, main="Before Normalization")
exprs_data <- normalizeBetweenArrays(exprs_data)
boxplot(exprs_data, outline = FALSE, las=2, main="After Normalization")

metadata <- pData(gse_data)
group <- ifelse(metadata$`endometriosis/no endometriosis:ch1` == "Endometriosis", "case", "control")
group <- factor(group)
table(group)

metadata <- pData(gse_data)

group <- ifelse(grepl("endometriosis", metadata$title, ignore.case = TRUE), "control", "endometriosis")
group <- factor(group)
table(group)

design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

fit <- lmFit(exprs_data, design)
contrast_matrix <- makeContrasts(endometriosis - control, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

results <- topTable(fit2, adjust = "fdr", number = Inf)
head(results)

gpl <- getGEO("GPL6480", AnnotGPL = TRUE)
annot_table <- Table(gpl)
annot <- annot_table[, c("ID", "Gene symbol")]
exprs_data$GeneSymbol <- annot$`Gene symbol`[match(rownames(exprs_data), annot$ID)]


deg <- subset(results, abs(logFC) > 1 & adj.P.Val < 0.05)
deg <- deg %>% filter(GeneSymbol != "")
head(deg)

write_xlsx(deg, "deg_GSE120103.xlsx")

deg_filtered <- deg %>% filter(!is.na(deg$GeneSymbol))

deg_filtered <- deg %>%
  filter(adj.P.Val < 0.05 & abs(logFC) > 1)
nrow(deg_filtered)

up_genes <- deg_filtered %>% filter(logFC> 1)
down_genes <- deg_filtered %>% filter(logFC < -1)

write_xlsx(up_genes, "deg_GSE120103_up_genes.xlsx")
write_xlsx(down_genes, "deg_GSE120103_down_genes.xlsx")

# Volkan Plotu
results$threshold <- as.factor(abs(results$logFC) > 1 & results$adj.P.Val < 0.05)

ggplot(results, aes(x = logFC, y = -log10(adj.P.Val), color = threshold)) +
  geom_point(alpha = 0.8) +
  theme_minimal() +
  scale_color_manual(values = c("grey", "red")) +
  labs(title = "Volcano Plot")

# Heatmap (ilk 50 DEG)
deg_rows <- rownames(deg_filtered)[1:50]
common_genes <- deg_rows[deg_rows %in% rownames(exprs_data)]
exprs_subset <- exprs_data[common_genes, ]
pheatmap(exprs_subset, scale = "row", main = "Heatmap of Top 50 DEGs")



deg_filtered$EntrezID <- mapIds(org.Hs.eg.db,
                                keys = deg_filtered$GeneSymbol,
                                column = "ENTREZID",
                                keytype = "SYMBOL",
                                multiVals = "first")

entrez_ids <- deg_filtered$EntrezID

#GO_BP
go_result_BP <- enrichGO(gene = entrez_ids,
                         OrgDb = org.Hs.eg.db,
                         keyType = "ENTREZID",
                         ont = "BP",
                         pAdjustMethod = "BH",
                         qvalueCutoff = 0.05,
                         readable = TRUE)

head(go_result_BP)
barplot(go_result_BP, showCategory = 20, title = "GO BP Enrichment")
write_xlsx(as.data.frame(go_result_BP), "deg_GSE120103_GO_BP_enrichment.xlsx")

#GO_MF
go_result_MF<- enrichGO(gene = entrez_ids,
                        OrgDb = org.Hs.eg.db,
                        keyType = "ENTREZID",
                        ont = "MF",
                        pAdjustMethod = "BH",
                        qvalueCutoff = 0.05,
                        readable = TRUE)

head(go_result_MF)
barplot(go_result_MF, showCategory = 20, title = "GO MF Enrichment")
write_xlsx(as.data.frame(go_result_MF), "deg_GSE120103_GO_MF_enrichment.xlsx")

#GO_CC
go_result_CC <- enrichGO(gene = entrez_ids,
                         OrgDb = org.Hs.eg.db,
                         keyType = "ENTREZID",
                         ont = "CC",
                         pAdjustMethod = "BH",
                         qvalueCutoff = 0.05,
                         readable = TRUE)

head(go_result_CC)
barplot(go_result_CC, showCategory = 20, title = "GO CC Enrichment")
write_xlsx(as.data.frame(go_result_CC), "deg_GSE120103_GO_CC_enrichment.xlsx")

#KEGG enrichment
kegg_enrich <- enrichKEGG(gene = entrez_ids,
                          organism = "hsa",
                          pAdjustMethod = "BH",
                          qvalueCutoff = 0.05)

kegg_enrich <- setReadable(kegg_enrich, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
head(kegg_enrich)
barplot(kegg_enrich, showCategory = 20, title = "KEGG Pathway Enrichment")

write_xlsx(as.data.frame(kegg_enrich), "deg_GSE120103_KEGG_enrichment.xlsx")

#Reactome Enrichment
reactome_enrich <- enrichPathway(gene = entrez_ids, 
                                 organism = "human", 
                                 pvalueCutoff = 0.1, 
                                 qvalueCutoff = 0.2)
barplot(reactome_enrich, showCategory = 20, title = "Reactome Pathway Enrichment")
write_xlsx(as.data.frame(reactome_enrich), "deg_GSE120103_reactome_enrichment.xlsx")





#xcell inf
exprs_data_symbols <- exprs_data
rownames(exprs_data_symbols) <- gene_symbols
exprs_data_symbols <- exprs_data_symbols[gene_symbols != "", ]


xcell_result <- xCellAnalysis(exprs_data_symbols)

xcell_df <- as.data.frame(xcell_result)
xcell_df <- cbind(Cell_Type = rownames(xcell_result), xcell_df)

pheatmap(xcell_result, 
         main = "xCell Immune Infiltration - GSE120103",
         scale = "row",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean")

write_xlsx(xcell_df, "xcell_GSE120103_results.xlsx")

