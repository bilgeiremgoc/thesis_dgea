library(GEOquery)
library(limma)
library(annotate)
library(ggplot2)
library(pheatmap)
library(hgu133a.db)
library(DESeq2)
library(edgeR)
library(writexl)

#Dataseti elde etme ve matrix oluşturma
gse <- getGEO("GSE7305")

gse1em <- getGEO("GSE7305", GSEMatrix = TRUE, AnnotGPL = TRUE)

exprs_data <- exprs(gse1em[[1]])

if (max(exprs_data, na.rm = TRUE) > 100) {
  exprs_data <- log2(exprs_data + 1)
}

#Normalizasyon var mı kontrol etme ve normalize etme
boxplot(exprs_data, outline = FALSE, las=2, main="Before Normalization")

exprs_data <- normalizeBetweenArrays(exprs_data)
boxplot(exprs_data, outline = FALSE, las=2, main="After Normalization")

#Control ve Hastalık gruplarını oluşturma
gse_data<- gse[[1]]

exprs(gse_data) <- exprs_data


metadata <- pData(gse_data)
group <- ifelse(grepl("normal", metadata$title, ignore.case = TRUE), "control", "endometriosis")
group <- factor(group)
table(group)  # grup sayısını kontrol et

design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)


fit <- lmFit(gse_data, design)
contrast_matrix <- makeContrasts(endometriosis - control, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

results <- topTable(fit2, adjust = "fdr", number = Inf)

# Volkan Plotu
results$threshold <- as.factor(abs(results$logFC) > 1 & results$adj.P.Val < 0.05)

ggplot(results, aes(x = logFC, y = -log10(adj.P.Val), color = threshold)) +
  geom_point(alpha = 0.8) +
  theme_minimal() +
  scale_color_manual(values = c("grey", "red")) +
  labs(title = "Volcano Plot")


deg <- subset(results, abs(logFC) > 1 & adj.P.Val < 0.05)
head(deg)
nrow(deg)


# Heatmap (ilk 50 DEG)
pheatmap(exprs_data[rownames(deg)[1:50], ], scale = "row")
#Annotation sorunu yaşandığı için çıkmadı ama gen listesi var


deg$GeneSymbol <- getSYMBOL(rownames(deg), "hgu133a.db")
head(deg)

deg_clean <- deg[!is.na(deg$GeneSymbol), ]


write_xlsx(deg_clean, "deg_GSE7305.xlsx")

deg_filtered <- deg %>%
  filter(adj.P.Val < 0.05 & abs(logFC) > 1)

# Kaç tane gen kaldı?
nrow(deg_filtered)

# İlk birkaç tanesini gör
head(deg_filtered)

up_genes <- deg_filtered %>% filter(logFC > 1)
down_genes <- deg_filtered %>% filter(logFC < -1)


