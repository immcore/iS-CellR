#!/usr/bin/env Rscript

if(isS4(scObject$val) || "defLabels" %in% isolate(input$clustLabels) || input$changeLabels)
{
	################## for Custom labels ################
  if("customLabels" %in% isolate(input$clustLabels)) {
    cluster.ids <- as.character(unlist(ClusterLabInfo$val[,1]))#, decreasing = FALSE)#c("CD4", "Bcells", "CD8cells",
    new.cluster.ids <- as.character(unlist(ClusterLabInfo$val[,2]))#c("CD4", "Bcells", "CD8cells", 
    if(input$changeLabels){
    	scObject$val@ident <- plyr::mapvalues(x = scObject$val@ident, from = cluster.ids, to = new.cluster.ids)
    }
    CInfo <- cbind(cluster.ids,new.cluster.ids)
    ClusterLabInfo$val <- CInfo
  }

  if("defLabels" %in% isolate(input$clustLabels)) {
    if(!is.null(dfcluster.ids$val)){
      new.cluster.ids <- as.character(unlist(ClusterLabInfo$val[,2]))
      current.ids <- as.character(unlist(ClusterLabInfo$val[,1]))#, decreasing = FALSE)#c("CD4", "Bcells", "CD8cells",
      scObject$val@ident <- plyr::mapvalues(x = scObject$val@ident, from = new.cluster.ids, to = current.ids)
    } else {
      new.cluster.ids = ""
      current.ids <- sort(as.character(unique(scObject$val@ident)), decreasing = FALSE)
    }
    cluster.ids <- current.ids
    CInfo <- cbind(cluster.ids,new.cluster.ids)
    ClusterLabInfo$val <- CInfo 
    dfcluster.ids$val <- cluster.ids      
  }

########### tSNE plot ggplot2 
# Create data frame of clusters computed by Seurat
df.cluster <- data.frame(Cell = names(scObject$val@ident), Cluster = scObject$val@ident)
    
# Create data frame of tSNE compute by Seurat
df.tsne <- data.frame(scObject$val@dr$umap@cell.embeddings)
# Add Cell column
df.tsne$Cell = rownames(df.tsne)
# Merge tSNE data frame to Cluster data frame
df.tsne <- merge(df.tsne, df.cluster, by = "Cell")

# Make df.tsne global 
tSNEmatrix$val <- df.tsne
mode$m <- 0
}

if("useheader" %in% isolate(input$clustLabels)) {
	current.clustID <- as.data.frame(scObject$val@ident)
	setDT(current.clustID, keep.rownames = TRUE)[]
	colnames(current.clustID) <- c("clust", "ident")
	current.clustID$clust <- gsub("\\..*|_.*|-.*", "", current.clustID$clust)
	current.clustID %>% group_by(ident,clust) %>% tally() -> lables
	new.lables <- as.data.frame(lables)
	new.lables <- new.lables[order(-new.lables$n), ]
	ClustLab <- subset(new.lables, !duplicated(ident))
	new.ident <- as.vector(ClustLab$ident)
	new.clust <- as.vector(ClustLab$clust)

	scObject$val@ident <- plyr::mapvalues(x = scObject$val@ident, from = new.ident, to = new.clust)
} 

#Five visualizations of marker gene expression
features.plot <- as.character(unlist(strsplit(input$DotplotGenes,","))) 
# each cluster

noGenes <- c()
Genes <- c()

for (i in features.plot) { 

	if(i %in% scObject$val@data@Dimnames[[1]])
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
	DownloadPlot$val$Dotplot <- DotPlot(object = scObject$val, genes.plot = Genes, plot.legend = TRUE, do.return = TRUE)
}



