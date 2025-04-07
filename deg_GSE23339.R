gse <- getGEO("GSE23339")

gse1em <- getGEO("GSE23339", GSEMatrix = TRUE, AnnotGPL = TRUE)

exprs_data <- exprs(gse1em[[1]])


# Veri hakkında genel bilgi
head(exprs(data))

boxplot(exprs_data, outline = FALSE, las=2, main="Before Normalization")


group <- factor(c(rep("endometriosis", 10), rep("control", 9)))
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
pheatmap(exprs_data[rownames(deg)[1:50], ], scale = "row")


annotation(gse1em[[1]])


feature_data <- fData(gse1em[[1]])
head(feature_data[, 1:5])  # İlk birkaç kolona bakalım, hangisi gen sembolü?

deg$GeneSymbol <- feature_data[rownames(deg), "Gene symbol"]
head(deg)

write.csv(deg, "deg_GSE23339.csv")

