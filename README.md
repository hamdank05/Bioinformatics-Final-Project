# Bioinformatics-Final-Project
## Background
###Code

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
library(rnaseqGene) # install using BiocManager
library(fission) # install using BiocManager
library(DESeq2)  # install using BiocManager
library(pheatmap)
library(RColorBrewer)
data(fission)
dds_o <- DESeqDataSet(fission, design = ~ strain + minute) #probably have to edit
rs <- rowSums(assay(fission))
smallestGroupSize <- 4
keep <- rowSums(counts(dds_o) >= 10) >= smallestGroupSize
dds <- dds_o[keep,]
means <- rowMeans(assay(dds))
sds <- rowSds(assay(dds))
hd <- tibble(m = means, s = sds)
vsd <- vst(dds, blind = FALSE)
meansv <- rowMeans(assay(vsd))
sdsv <- rowSds(assay(vsd))
hdv <- data_frame(m = meansv, s = sdsv)
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( vsd$minute, vsd$strain, sep = " - " )
colnames(sampleDistMatrix) <- paste()
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors, show_colnames = T, show_rownames = T)
