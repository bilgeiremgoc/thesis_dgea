# Gerekli kütüphaneler
library(ggplot2)
library(pheatmap)
library(dplyr)

gse <- getGEO("GSE7307")

gse1em <- getGEO("GSE7307", GSEMatrix = TRUE, AnnotGPL = TRUE)

exprs_data <- exprs(gse1em[[1]])


if (max(exprs_data, na.rm = TRUE) > 100) {
  exprs_data <- log2(exprs_data + 1)
}


boxplot(exprs_data, outline = FALSE, las=2, main="Before Normalization")

exprs_data <- normalizeBetweenArrays(exprs_data)
boxplot(exprs_data, outline = FALSE, las=2, main="After Normalization")

gse_data<- gse[[1]]

metadata <- pData(gse_data)


head(metadata$characteristics_ch1)

# Endometrium doku olanları seç
endometrium_samples <- grepl("endometrium", metadata$characteristics_ch1, ignore.case = TRUE)

# Veriyi filtrele
exprs_data <- exprs_data[, endometrium_samples]
metadata <- metadata[endometrium_samples, ]

group <- ifelse(grepl("normal", metadata$`Disease type:ch1`, ignore.case = TRUE), "endometriosis", "normal")
group <- factor(group)
table(group)

design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

fit <- lmFit(exprs_data, design)
contrast_matrix <- makeContrasts(endometriosis - normal, levels = design)
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

deg <- deg[!is.na(deg$GeneSymbol) & deg$GeneSymbol != "", ]
head(deg)

write_xlsx(deg, "deg_GSE7307.xlsx")

deg_filtered <- deg %>%
  filter(adj.P.Val < 0.05 & abs(logFC) > 1)

# Kaç tane gen kaldı?
nrow(deg_filtered)

# İlk birkaç tanesini gör
head(deg_filtered)

up_genes <- deg_filtered %>% filter(logFC > 1)
down_genes <- deg_filtered %>% filter(logFC < -1)
