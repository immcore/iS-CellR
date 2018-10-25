#!/usr/bin/env Rscript 

if(input$SwitchUpload == "TRUE"){
  inFile <- input$file1
} else {
  inFile <- parseFilePaths(roots=c(root='.'), input$sfile)
}

if(file_ext(inFile$name) == "rds"){

	scObject$val <- countRds$val
	countFileQC$val <- as.data.frame(as.matrix(scObject$val@data))
	
} else {

	countFile <- data.matrix(countFile$val)
   	
   	scObject$val <- CreateSeuratObject(raw.data = countFile, min.cells = 3, min.genes = 3, 
                                  project = values$ProjectName)

	countFileQC$val <- as.data.frame(as.matrix(scObject$val@data))

	# Normalizing the data
	scObject$val <- NormalizeData(object = scObject$val, normalization.method = "LogNormalize", 
	                             scale.factor = 10000)
	    
	# Detection of variable genes across the single cells
	scObject$val <- FindVariableGenes(object = scObject$val, mean.function = ExpMean, dispersion.function = LogVMR, 
		x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, do.plot = FALSE)

	hv.genes <- head(rownames(scObject$val@hvg.info), 1000)

	# Scaling the data and removing unwanted sources of variation
	mito.genes <- grep(pattern = "^mt-", ignore.case = TRUE, x = rownames(x = scObject$val@data), value = TRUE)
	percent.mito <- Matrix::colSums(scObject$val@raw.data[mito.genes, ])/Matrix::colSums(scObject$val@raw.data)
	scObject$val <- AddMetaData(object = scObject$val, metadata = percent.mito, col.name = "percent.mito")

	ribo.genes <- grep(pattern = "^rps|^rpl", ignore.case = TRUE, x = rownames(x = scObject$val@data), value = TRUE)
	percent.ribo <- Matrix::colSums(scObject$val@raw.data[ribo.genes, ])/Matrix::colSums(scObject$val@raw.data)
	scObject$val <- AddMetaData(object = scObject$val, metadata = percent.ribo, col.name = "percent.ribo")

	scObject$val <- ScaleData(object = scObject$val, genes.use = hv.genes, display.progress = FALSE, 
	    vars.to.regress = c("percent.mito","percent.ribo"), do.par = TRUE, num.cores = 2)

	# Perform linear dimensional reduction
	scObject$val <- RunPCA(object = scObject$val, pc.genes = hv.genes, pcs.compute = 100, do.print = FALSE, 
	    pcs.print = 1:5, genes.print = 5)

	scObject$val <- FindClusters(object = scObject$val, reduction.type = "pca", dims.use = 1:75, resolution = 0.9, 
	    save.SNN = FALSE, n.start = 10, nn.eps = 0.5, print.output = FALSE)

	scObject$val <- RunTSNE(object = scObject$val, reduction.use = "pca", dims.use = 1:75, tsne.method = "FIt-SNE", nthreads = 4, reduction.name = "FItSNE", reduction.key = "FItSNE_", fast_tsne_path = "/bin/fast_tsne", max_iter = 2000)
	scObject$val <- RunUMAP(object = scObject$val, reduction.use = "pca", dims.use = 1:75, min_dist = 0.75)
} 

	mode$n <- 0
	mode$m <- 1
	mode$l <- 1
