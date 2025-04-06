library(GEOquery)
library(limma)
library(annotate)
library(ggplot2)
library(pheatmap)
library(hgu133a.db)

#Datasetin database üzerinden çekilmesi

gse <- getGEO("GSE25628", GSEMatrix = TRUE)
exprSet <- exprs(gse[[1]])
phenoData <- pData(gse[[1]])

#Data normalize edilmiş mi kontrol edilmesi

boxplot(exprSet, outline = FALSE, las=2, main="Before Normalization")

exprSet <- normalizeBetweenArrays(exprSet)
boxplot(exprSet, outline = FALSE, las=2, main="After Normalization")


