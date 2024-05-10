# Bioinformatics-Final-Project
## Background

Understanding how organisms respond to environmental stresses at the molecular level is fundamental to unraveling cellular mechanisms and adaptation strategies. Here, we investigate the gene expression dynamics of Schizosaccharomyces pombe (fission yeast) in response to oxidative stress.

Two strains of yeast were exposed to 1M sorbitol treatment for different time periods. The strains consisted of an unedited wild type and a mutant with deletion of atf21, a key transcription factor implicated in stress response pathways. The strains were exposed to treatment for 0, 15, 30, 60, 120, and 180 minutes.

Here we show the correlation of gene expression under different time periods of oxidative stress from the wt and mutant strains through the creation of a heatmap. 


### Code

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
  
library(rnaseqGene) # install using BiocManager

library(fission) # install using BiocManager

library(DESeq2)  # install using BiocManager

library(pheatmap)

library(RColorBrewer)

data(fission)

dds_o <- DESeqDataSet(fission, design = ~ strain + minute) 

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
