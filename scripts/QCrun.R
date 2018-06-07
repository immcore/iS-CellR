#!/usr/bin/env Rscript

# select any cells and rows with expr > 0
countFile <- data.matrix(countFile$val[apply(countFile$val[, -1], MARGIN = 1, function(x) any(x > input$Expr2)), ])

countMatrix <- CreateSeuratObject(raw.data = countFile, min.cells = input$nCells, min.genes = input$nGenes2, 
                                  project = values$ProjectName)
    
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