# Gerekli kütüphaneler
library(ggplot2)
library(pheatmap)
library(dplyr)

deg_data <- read.table("GSE7307.top.table.tsv", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)

nrow(deg_data)

# İlk satırlara göz at
head(deg_data)

# Örnek eşik: logFC > 1 veya < -1 ve adj.P.Val < 0.05
deg_filtered <- deg_data %>%
  filter(adj.P.Val < 0.05 & abs(logFC) > 1)

# Kaç tane gen kaldı?
nrow(deg_filtered)

# İlk birkaç tanesini gör
head(deg_filtered)

up_genes <- deg_filtered %>% filter(logFC > 1)
down_genes <- deg_filtered %>% filter(logFC < -1)

cat("Upregulated:", nrow(up_genes), "\nDownregulated:", nrow(down_genes), "\n")
