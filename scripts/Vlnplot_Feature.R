#!/usr/bin/env Rscript

if(!isS4(tSNEClusters$val) || (mode$n == 0 && mode$m == 1) || "defLabels" %in% isolate(input$clustLabels) || input$changeLabels)
{
	if(!isS4(tSNEClusters$val) || (mode$n == 0 && mode$m == 1))
	{
		Vlndata <- ProjectPCA(object = seuratObject$val, do.print = FALSE)
 		Vlndata <- FindClusters(object = Vlndata, reduction.type = "pca", dims.use = 1:10, 
                            resolution = 0.6, print.output = 0, save.SNN = TRUE)
		# Make tsne global 
		tSNEClusters$val <- Vlndata
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

#Five visualizations of marker gene expression
features.plot <- as.character(unlist(strsplit(input$VlnplotGenes," "))) 
# Violin plots. Visualize single cell expression distributions in each
# cluster

########### tSNE plot ggplot2 
# Create data frame of clusters computed by Seurat

plot_list <- list()
noGenes <- c()

for (i in features.plot) { 

	if(i %in% tSNEClusters$val@data@Dimnames[[1]])
	{
		df <- data.frame(Cell = names(tSNEClusters$val@ident), Expression = tSNEClusters$val@data[i,], Gene = i)
		df.cluster <- data.frame(Cell = names(tSNEClusters$val@ident), Cluster = tSNEClusters$val@ident)
		df.vln <- merge(df, df.cluster, by = "Cell")

		if("useheader" %in% isolate(input$clustLabels)) {
			df.vln$Celltype <- df.vln$Cell
			df.vln$Celltype <- gsub("\\..*|_.*|-.*", "", df.vln$Celltype)

			plot_single <- ggplot(data = df.vln, aes(x=factor(x=Celltype), y=Expression)) + geom_violin(scale = "width", adjust.use = 1, trim=TRUE, mapping = aes(fill = factor(x = Celltype))) + 
					geom_jitter(height = 0, size = 1) + xlab("Cluster") + ylab("Expression") + 
					theme(axis.text.x = element_text(angle = 35, vjust=1.0, hjust=1.0)) +
					theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
					ggtitle(i) + theme(plot.title = element_text(color="black", size=14, face="bold")) + theme(legend.position="none") + theme(legend.title = element_blank())
		} else {
			plot_single <- ggplot(data = df.vln, aes(x=factor(x=Cluster), y=Expression)) + geom_violin(scale = "width", adjust.use = 1, trim=TRUE, mapping = aes(fill = factor(x = Cluster))) + 
					geom_jitter(height = 0, size = 1) + xlab("Cluster") + ylab("Expression") + 
					theme(axis.text.x = element_text(angle = 35, vjust=1.0, hjust=1.0)) +
					theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
					ggtitle(i) + theme(plot.title = element_text(color="black", size=14, face="bold")) + theme(legend.position="none") + theme(legend.title = element_blank())
		}

		plot_list[[length(plot_list)+1]] <- plot_single
	}
	else{
		noGenes[length(noGenes)+1] = i
	}
}

GenesAbsent$val <- noGenes

if(length(plot_list) == 1){
	num_col = 1
} else if(length(plot_list) == 2){
	num_col = 2
} else if(length(plot_list) >= 3 && length(plot_list) <= 9){
	num_col = 3
} else if(length(plot_list) > 9){
	num_col = 4
}

if(length(GenesAbsent$val) == length(features.plot)) {
      plot.new()
      DownloadPlot$val$Vlnplot <- NULL
      } else {
      	grid.arrange(grobs=plot_list, ncol=num_col)
      	DownloadPlot$val$Vlnplot <- plot_grid(plotlist = plot_list, ncol = num_col)
      }
