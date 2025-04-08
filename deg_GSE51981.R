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


# Veri hakkında genel bilgi
head(exprs_data)

boxplot(exprs_data, outline = FALSE, las=2, main="Before Normalization")

exprs_data <- normalizeBetweenArrays(exprs_data, method = "quantile")
boxplot(exprs_data, outline = FALSE, las=2, main="After Normalization")


group <- factor(c(rep("endometriosis", 77), rep("control", 71)))
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group) 


fit <- lmFit(exprs_data, design)

# Grup kar????la??t??rmas??
contrast_matrix <- makeContrasts(endometriosis - control, levels = design)
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

# Heatmap (ilk 50 DEG)
pheatmap(exprs_data[rownames(deg)[1:50], ], scale = "row")

deg$GeneSymbol <- getSYMBOL(rownames(deg), "hgu133a.db")
head(deg)

write.csv(deg, "deg_GSE51981.csv")


deg_filtered <- deg %>%
  filter(adj.P.Val < 0.05 & abs(logFC) > 1)

up_genes <- deg_filtered %>% filter(logFC > 1)
down_genes <- deg_filtered %>% filter(logFC < -1)

cat("Upregulated:", nrow(up_genes), "\nDownregulated:", nrow(down_genes), "\n")
