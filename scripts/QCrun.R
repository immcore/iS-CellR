#!/usr/bin/env Rscript

# select any cells and rows with expr > 0
countFile <- data.matrix(countFile$val[apply(countFile$val[, -1], MARGIN = 1, function(x) any(x > input$Expr2)), ])

scObject$val <- CreateSeuratObject(counts = countFile, min.cells = input$nCells, min.features = input$nGenes2, project = values$ProjectName)

countFileQC$val <- as.data.frame(as.matrix(scObject$val@assays$RNA@data), 'sparseMatrix')
         
# Scaling the data and removing unwanted sources of variation
mito.genes <- grep(pattern = "^mt-", ignore.case = TRUE, x = rownames(x = scObject$val), value = TRUE)
percent.mito <- Matrix::colSums(x = GetAssayData(object = scObject$val, slot = "counts")[mito.genes, ]) / Matrix::colSums(x = GetAssayData(object = scObject$val, slot = "counts"))
scObject$val[['percent.mito']] <- percent.mito
   
ribo.genes <- grep(pattern = "^rps|^rpl", ignore.case = TRUE, x = rownames(x = scObject$val), value = TRUE)
percent.ribo <- Matrix::colSums(x = GetAssayData(object = scObject$val, slot = "counts")[ribo.genes, ]) / Matrix::colSums(x = GetAssayData(object = scObject$val, slot = "counts"))
scObject$val[['percent.ribo']] <- percent.ribo
      
# Normalizing the data
scObject$val <- NormalizeData(object = scObject$val, normalization.method = "LogNormalize", scale.factor = 10000)
     
# Detection of variable genes across the single cells
scObject$val <- FindVariableFeatures(object = scObject$val, selection.method = 'mean.var.plot', mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(0.05, Inf))
          
scObject$val <- ScaleData(object = scObject$val, features = rownames(x = scObject$val), do.par = TRUE, num.cores = 4)#, vars.to.regress = c("percent.mito", "percent.ribo"), do.par = TRUE, num.cores = 4)
     
# Perform linear dimensional reduction
scObject$val <- RunPCA(object = scObject$val)
    
scObject$val <- ProjectDim(object = scObject$val)
    
#scObject$val <- JackStraw(object = scObject$val, num.replicate = 100)
#scObject$val <- ScoreJackStraw(object = scObject$val, dims = 1:20)
      
scObject$val <- FindNeighbors(object = scObject$val, dims = 1:10)
scObject$val <- FindClusters(object = scObject$val, resolution = 0.6)
      
scObject$val <- RunTSNE(object = scObject$val, dims = 1:10)
scObject$val <- RunTSNE(object = scObject$val, reduction.use = "pca", dims.use = 1:75, tsne.method = "FIt-SNE", nthreads = 4, reduction.name = "FItSNE", reduction.key = "FItSNE_", fast_tsne_path = "FIt-SNE-master/bin/fast_tsne", max_iter = 2000)
scObject$val <- RunUMAP(object = scObject$val, dims = 1:10)

mode$n <- 0
mode$m <- 1
mode$l <- 1
