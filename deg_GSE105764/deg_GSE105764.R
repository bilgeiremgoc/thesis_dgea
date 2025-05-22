gse <- getGEO("GSE105764", GSEMatrix = TRUE)
pheno <- pData(gse[[1]])

url <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE105nnn/GSE105764/suppl/GSE105764_count.txt.gz"
download.file(url, destfile = "GSE105764_count.txt.gz")

counts <- read.table("GSE105764_count.txt.gz", header = TRUE, row.names = 1)

print(dim(counts)) 
print(head(counts)) 

group <- ifelse(grepl("case", pheno$characteristics_ch1, ignore.case = TRUE), "control", "case")
group <- factor(group)
table(group)  


dge <- DGEList(counts = counts, group = group)

keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes=FALSE]

dge <- calcNormFactors(dge)

design <- model.matrix(~group)
dge <- estimateDisp(dge, design)

fit <- glmFit(dge, design)
lrt <- glmLRT(fit, coef=2)

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

deg_filtered <- deg_results[deg_results$FDR < 0.05 & abs(deg_results$logFC) > 1, ]
head(deg_filtered)

cat("AnlamlD1 gen sayD1sD1:", nrow(deg_filtered), "\n")

write_xlsx(deg_filtered, "deg_GSE105764.xlsx")

up_genes <- deg_filtered %>% filter(logFC > 1)
down_genes <- deg_filtered %>% filter(logFC < -1)

write_xlsx(up_genes, "deg_GSE105764_up_genes.xlsx")
write_xlsx(down_genes, "deg_GSE105764_down_genes.xlsx")


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

#Heatmap
top_genes <- rownames(deg_filtered)[1:50]
heatmap_data <- counts[top_genes, ]
heatmap_data <- log2(heatmap_data + 1)

pheatmap(heatmap_data, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         scale = "row", 
         show_rownames = TRUE,
         main = "Top 20 DEGs")


gene_list <- rownames(deg_results)
gene_df <- bitr(gene_list, fromType="ENSEMBL", toType=c("SYMBOL", "ENTREZID"), OrgDb="org.Hs.eg.db")
head(gene_df)



#GO_BP
go_result_BP <- enrichGO(gene = gene_df$ENTREZID,
                         OrgDb = org.Hs.eg.db,
                         ont = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.05,
                         readable = TRUE)

head(go_result_BP)
barplot(go_result_BP, showCategory = 20, title = "GO BP Enrichment")
write_xlsx(as.data.frame(go_result_BP), "deg_GSE105764_GO_BP_enrichment.xlsx")

#GO_MF
go_result_MF<- enrichGO(gene = gene_df$ENTREZID,
                        OrgDb = org.Hs.eg.db,
                        keyType = "ENTREZID",
                        ont = "MF",
                        pAdjustMethod = "BH",
                        qvalueCutoff = 0.05,
                        readable = TRUE)

head(go_result_MF)
barplot(go_result_MF, showCategory = 20, title = "GO MF Enrichment")
write_xlsx(as.data.frame(go_result_MF), "deg_GSE105764_GO_MF_enrichment.xlsx")

#GO_CC
go_result_CC <- enrichGO(gene = gene_df$ENTREZID,
                         OrgDb = org.Hs.eg.db,
                         keyType = "ENTREZID",
                         ont = "CC",
                         pAdjustMethod = "BH",
                         qvalueCutoff = 0.05,
                         readable = TRUE)

head(go_result_CC)
barplot(go_result_CC, showCategory = 20, title = "GO CC Enrichment")
write_xlsx(as.data.frame(go_result_CC), "deg_GSE105764_GO_CC_enrichment.xlsx")

#KEGG enrichment
kegg_enrich <- enrichKEGG(gene = gene_df$ENTREZID,
                          organism = "hsa",
                          pAdjustMethod = "BH",
                          qvalueCutoff = 0.05)

kegg_enrich <- setReadable(kegg_enrich, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
head(kegg_enrich)
barplot(kegg_enrich, showCategory = 20, title = "KEGG Pathway Enrichment")

write_xlsx(as.data.frame(kegg_enrich), "deg_GSE105764_KEGG_enrichment.xlsx")

#Reactome Enrichment
reactome_enrich <- enrichPathway(gene = gene_df$ENTREZID, 
                                 organism = "human", 
                                 pvalueCutoff = 0.1, 
                                 qvalueCutoff = 0.2)
barplot(reactome_enrich, showCategory = 20, title = "Reactome Pathway Enrichment")

write_xlsx(as.data.frame(reactome_enrich), "deg_GSE105764_reactome_enrichment.xlsx")


exprs_data <- cpm(dge, log = TRUE, prior.count = 1) 
gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = rownames(exprs_data),
                       column = "SYMBOL",
                       keytype = "ENSEMBL",
                       multiVals = "first")

exprs_data_symbols <- exprs_data
rownames(exprs_data_symbols) <- gene_symbols
exprs_data_symbols <- exprs_data_symbols[!is.na(rownames(exprs_data_symbols)), ]

exprs_data_symbols <- rowsum(exprs_data_symbols, group = rownames(exprs_data_symbols))


xcell_result <- xCellAnalysis(exprs_data_symbols)

xcell_df <- as.data.frame(xcell_result)
xcell_df <- cbind(Cell_Type = rownames(xcell_result), xcell_df)

# Is?? haritas??
pheatmap(xcell_result,
         main = "xCell Immune Infiltration - GSE105764",
         scale = "row",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean")


write_xlsx(xcell_df, "deg_GSE105764_xCell_Results.xlsx")



xcell_matrix <- as.data.frame(t(xcell_result))
xcell_matrix$sample <- rownames(xcell_matrix)
xcell_matrix$group <- group

long_xcell <- xcell_matrix %>%
  pivot_longer(cols = -c(sample, group), names_to = "cell_type", values_to = "proportion")

# Wilcoxon testleri
stat_tests <- long_xcell %>%
  group_by(cell_type) %>%
  summarise(
    p_value = wilcox.test(proportion ~ group)$p.value,
    mean_case = mean(proportion[group == "case"]),
    mean_control = mean(proportion[group == "control"])
  ) %>%
  mutate(adj_p = p.adjust(p_value, method = "fdr")) %>%
  arrange(p_value)

write_xlsx(stat_tests, "deg_GSE105764_xcell_stat_tests.xlsx")

ggplot(long_xcell, aes(x = group, y = proportion, fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~ cell_type, scales = "free_y") +
  theme_minimal() +
  labs(title = "Immune Cell Infiltration (xCell) - GSE105764") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

