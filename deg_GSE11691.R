install.packages("BiocManager")
BiocManager::install("GEOquery")
library(GEOquery)

BiocManager::install(c("GEOquery", "limma", "annotate", "hgu133a.db"))
install.packages("ggplot2")
install.packages("pheatmap")

library(GEOquery)
library(limma)
library(annotate)
library(ggplot2)
library(pheatmap)
library(hgu133a.db)

getGEO("GSE11691")
gse1em <- getGEO("GSE11691", GSEMatrix = TRUE, AnnotGPL = TRUE)

exprs_data <- exprs(gse1em[[1]])
exprs_data <- normalizeBetweenArrays(exprs_data, method = "quantile")

design <- model.matrix(~ 0 + factor(c(rep("endo", 9), rep("control", 9))))  
fit <- lmFit(exprs_data, design)

group <- factor(c(rep("endometriosis", 9), rep("control", 9)))
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)


boxplot(exprs_data, outline = FALSE, las=2, main="Before Normalization")

exprs_data <- normalizeBetweenArrays(exprs_data)
boxplot(exprs_data, outline = FALSE, las=2, main="After Normalization")


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

write.csv(deg, "deg_GSE11691.csv")








