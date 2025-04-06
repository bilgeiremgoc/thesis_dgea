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

gse <- getGEO("GSE25628", GSEMatrix = TRUE)

exprSet <- exprs(gse[[1]])  # Expression matrix
metaData <- pData(gse[[1]])

head(exprSet)
head(metaData)


# Grup bilgilerini düzenleyelim
metaData$disease <- factor(metaData$disease, levels = c("Normal (control)", "Normal (eutopic)", "Pathological (ectopic)"))

# Geçerli isimler oluşturmak için boşluk ve özel karakterleri değiştirelim
metaData$disease <- gsub(" ", ".", metaData$disease)  # Boşlukları nokta ile değiştirme
metaData$disease <- make.names(metaData$disease)  # Geçerli R isimleri yapmak

# Yeni isimlere göz atalım
table(metaData$disease)


group <- factor(c(rep("endometriosis", 7), rep("control", 6), rep("normal", 9)))

exprSet <- as.matrix(exprSet)

design <- model.matrix(~ 0 + group)
colnames(design) <- levels(metaData$disease)

# Normalizasyon için voom kullanıyoruz
y <- DGEList(counts = exprSet)
y <- calcNormFactors(y)
v <- voom(y, design = model.matrix(~disease, data = metaData), plot = TRUE)

# Linear model kurun
fit <- lmFit(v, design = model.matrix(~disease, data = metaData))

fit <- lmFit(exprSet, design)

contrast.matrix <- makeContrasts(
  Normal.control - Normal.eutopic,       
  Normal.control - Pathological.ectopic, 
  Normal.eutopic - Pathological.ectopic, 
  levels = design
)

fit2 <- contrasts.fit(fit, contrast.matrix)

# Sonuçları değerlendirin
fit2 <- eBayes(fit2)
