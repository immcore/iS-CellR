#!/usr/bin/env Rscript

countFile <- data.matrix(countFile$val)
   
countMatrix <- CreateSeuratObject(raw.data = countFile, min.cells = 3, min.genes = 3, 
                                  project = "")

countFileQC$val <- as.data.frame(countMatrix@data)

# Normalizing the data
countMatrix <- NormalizeData(object = countMatrix, normalization.method = "LogNormalize", 
                             scale.factor = 10000)
    
# Detection of variable genes across the single cells
countMatrix <- FindVariableGenes(object = countMatrix, mean.function = ExpMean, dispersion.function = LogVMR, 
                                 x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, do.plot = FALSE)
   
countMatrix <- ScaleData(object = countMatrix)
   
# Perform linear dimensional reduction
seuratObject$val <- RunPCA(object = countMatrix, pc.genes = countMatrix@var.genes, do.print = TRUE, pcs.print = 1:5, 
                      genes.print = 5)

mode$n <- 0
mode$m <- 1
mode$l <- 1