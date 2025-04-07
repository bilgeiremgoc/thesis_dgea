library(GEOquery)
library(limma)
library(annotate)
library(ggplot2)
library(pheatmap)
library(hgu133a.db)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
library(DESeq2)

BiocManager::install("edgeR")
library(edgeR)

gse <- getGEO("GSE25628")

gse1em <- getGEO("GSE2628", GSEMatrix = TRUE, AnnotGPL = TRUE)

exprs_data <- exprs(gse1em[[1]])

head(exprs(data))

boxplot(exprs_data, outline = FALSE, las=2, main="Before Normalization")

exprs_data <- normalizeBetweenArrays(exprs_data)
boxplot(exprs_data, outline = FALSE, las=2, main="After Normalization")

group <- factor(c(rep("endometriosis", 7), rep("control", 15)))
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

ggplot(results, aes(x = logFC, y = -log10(adj.P.Val), color = threshold)) +
  geom_point(alpha = 0.8) +
  theme_minimal() +
  scale_color_manual(values = c("grey", "red")) +
  labs(title = "Volcano Plot")


nrow(deg)


# Heatmap (ilk 50 DEG)
pheatmap(exprs_data[rownames(deg)[1:50], ], scale = "row")
#Annotation sorunu yaşandığı için çıkmadı ama gen listesi var


deg$GeneSymbol <- getSYMBOL(rownames(deg), "hgu133a.db")
head(deg)

write.csv(deg, "deg_GSE25628.csv")

head(rownames(exprs_data))
head(rownames(deg))


