#!/usr/bin/env Rscript

if(!isS4(tSNEClusters$val) || (mode$n == 0 && mode$m == 1) || "defLabels" %in% isolate(input$clustLabels) || input$changeLabels)
{
		if(!isS4(tSNEClusters$val) || (mode$n == 0 && mode$m == 1))
	{
		Dotdata <- ProjectPCA(object = seuratObject$val, do.print = FALSE)
 		Dotdata <- FindClusters(object = Dotdata, reduction.type = "pca", dims.use = 1:10, 
                            resolution = 0.6, print.output = 0, save.SNN = TRUE)
 		# Make tsne global 
		#assign('tSNEClusters', Dotdata, envir=.GlobalEnv)
		tSNEClusters$val <- Dotdata
	}
	
	################## for Custom labels ################
  if("customLabels" %in% isolate(input$clustLabels)) {
    cluster.ids <- as.character(unlist(ClusterLabInfo$val[,1]))#, decreasing = FALSE)#c("CD4", "Bcells", "CD8cells",
    new.cluster.ids <- as.character(unlist(ClusterLabInfo$val[,2]))#c("CD4", "Bcells", "CD8cells", 
    tSNEClusters$val@ident <- plyr::mapvalues(x = tSNEClusters$val@ident, from = cluster.ids, to = new.cluster.ids)
    CInfo <- cbind(cluster.ids,new.cluster.ids)
    ClusterLabInfo$val <- CInfo
  }

  if("defLabels" %in% isolate(input$clustLabels)) {
    if(!is.null(dfcluster.ids$val)){
      new.cluster.ids <- as.character(unlist(ClusterLabInfo$val[,2]))
      current.ids <- as.character(unlist(ClusterLabInfo$val[,1]))#, decreasing = FALSE)#c("CD4", "Bcells", "CD8cells",
      tSNEClusters$val@ident <- plyr::mapvalues(x = tSNEClusters$val@ident, from = new.cluster.ids, to = current.ids)
    } else {
      new.cluster.ids = ""
      current.ids <- sort(as.character(unique(tSNEClusters$val@ident)), decreasing = FALSE)
    }
    cluster.ids <- current.ids
    CInfo <- cbind(cluster.ids,new.cluster.ids)
    ClusterLabInfo$val <- CInfo 
    dfcluster.ids$val <- cluster.ids      
  }

mode$m <- 0
}

if("useheader" %in% isolate(input$clustLabels)) {
	current.clustID <- as.data.frame(tSNEClusters$val@ident)
	setDT(current.clustID, keep.rownames = TRUE)[]
	colnames(current.clustID) <- c("clust", "ident")
	current.clustID$clust <- gsub("\\..*|_.*|-.*", "", current.clustID$clust)
	current.clustID %>% group_by(ident,clust) %>% tally() -> lables
	new.lables <- as.data.frame(lables)
	new.lables <- new.lables[order(-new.lables$n), ]
	ClustLab <- subset(new.lables, !duplicated(ident))
	new.ident <- as.vector(ClustLab$ident)
	new.clust <- as.vector(ClustLab$clust)

	tSNEClusters$val@ident <- plyr::mapvalues(x = tSNEClusters$val@ident, from = new.ident, to = new.clust)
} 

#Five visualizations of marker gene expression
features.plot <- as.character(unlist(strsplit(input$DotplotGenes,","))) 

# Joy plots - from ggjoy. Visualize single cell expression distributions in
# each cluster

noGenes <- c()
Genes <- c()

for (i in features.plot) { 

	if(i %in% tSNEClusters$val@data@Dimnames[[1]])
	{
		Genes[length(Genes)+1] = i
	}
	else{
		noGenes[length(noGenes)+1] = i
	}
}

#assign('GenesAbsent', noGenes, envir=.GlobalEnv)
GenesAbsent$val <- noGenes

if(length(GenesAbsent$val) == length(features.plot)) {
	plot.new()
	DownloadPlot$val$Dotplot <- NULL
} else {
	DownloadPlot$val$Dotplot <- DotPlot(object = tSNEClusters$val, genes.plot = Genes, plot.legend = TRUE, do.return = TRUE)
}



