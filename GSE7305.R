gse <- getGEO("GSE7305")

gse1em <- getGEO("GSE7305", GSEMatrix = TRUE, AnnotGPL = TRUE)

exprs_data <- exprs(gse1em[[1]])


# Veri hakkında genel bilgi
head(exprs(data))


boxplot(exprs_data, outline = FALSE, las=2, main="Before Normalization")

exprs_data <- normalizeBetweenArrays(exprs_data)
boxplot(exprs_data, outline = FALSE, las=2, main="After Normalization")

group <- factor(c(rep("endometriosis", 10), rep("control", 10)))
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group) 


fit <- lmFit(exprs_data, design)
contrast_matrix <- makeContrasts(endometriosis - control, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)


results <- topTable(fit2, adjust = "fdr", number = Inf)

deg <- subset(results, abs(logFC) > 1 & adj.P.Val < 0.05)
head(deg)

# Volkan Plotu
results$threshold <- as.factor(abs(results$logFC) > 1 & results$adj.P.Val < 0.05)



# Heatmap (ilk 50 DEG)
pheatmap(exprs_data[rownames(deg)[1:50], ], scale = "row")

deg$GeneSymbol <- getSYMBOL(rownames(deg), "hgu133a.db")
head(deg)

deg$Significant <- "Not Significant"
deg$Significant[deg$adj.P.Val < 0.05 & abs(deg$logFC) > 0.5] <- "Significant"

head(deg)  # İlk birkaç satırı gözden geçir
colnames(deg)  # Sütun isimlerini kontrol et


ggplot(deg, aes(x = logFC, y = -log10(adj.P.Val), color = Significant)) +
  geom_point(alpha = 0.5, size = 2) +
  scale_color_manual(values = c("Not Significant" = "gray", "Significant" = "red")) +
  theme_minimal() +
  theme(legend.position = "top") +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 Adjusted p-value")


write.csv(deg, "deg_GSE7305.csv")


