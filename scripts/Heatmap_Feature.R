#!/usr/bin/env Rscript

if(isS4(scObject$val) || "defLabels" %in% isolate(input$clustLabels) || input$changeLabels)
{
	################## for Custom labels ################
  if("customLabels" %in% isolate(input$clustLabels)) {
    cluster.ids <- as.character(unlist(ClusterLabInfo$val[,1]))#, decreasing = FALSE)#c("CD4", "Bcells", "CD8cells",
    new.cluster.ids <- as.character(unlist(ClusterLabInfo$val[,2]))#c("CD4", "Bcells", "CD8cells", 
    if(input$changeLabels){
      #scObject$val <- RenameIdents(scObject$val, new.cluster.ids)
      Idents(scObject$val) <- plyr::mapvalues(x = Idents(scObject$val), from = cluster.ids, to = new.cluster.ids)
    }
    CInfo <- cbind(cluster.ids,new.cluster.ids)
    ClusterLabInfo$val <- CInfo
  }

  if("defLabels" %in% isolate(input$clustLabels)) {
    if(!is.null(dfcluster.ids$val)){
      new.cluster.ids <- as.character(unlist(ClusterLabInfo$val[,2]))
      current.ids <- as.character(unlist(ClusterLabInfo$val[,1]))#, decreasing = FALSE)#c("CD4", "Bcells", "CD8cells",
      #scObject$val <- RenameIdents(scObject$val, current.ids)
      Idents(scObject$val) <- plyr::mapvalues(x = Idents(scObject$val), from = new.cluster.ids, to = current.ids)
    } else {
      new.cluster.ids = ""
      current.ids <- sort(as.character(unique(Idents(scObject$val))), decreasing = FALSE)
    }
    cluster.ids <- current.ids
    CInfo <- cbind(cluster.ids,new.cluster.ids)
    ClusterLabInfo$val <- CInfo 
    dfcluster.ids$val <- cluster.ids      
  }

########### tSNE plot ggplot2 
# Create data frame of clusters computed by Seurat
df.cluster <- data.frame(Cell = names(Idents(object = scObject$val)), Cluster = Idents(object = scObject$val))
    
# Create data frame of tSNE compute by Seurat
df.umap <- data.frame(Embeddings(object = scObject$val, reduction = "umap"))
# Add Cell column
colnames(df.umap) <- c("UMAP1","UMAP2")
df.umap$Cell = rownames(df.umap)
# Create data frame of tSNE compute by Seurat
df.FItsne <- data.frame(Embeddings(object = scObject$val, reduction = "FItSNE"))
# Add Cell column
df.FItsne$Cell = rownames(df.FItsne)

# Merge tSNE data frame to Cluster data frame
df.tsne <- merge(df.umap, df.FItsne, by = "Cell")
df.tsne <- merge(df.tsne, df.cluster, by = "Cell")

# Make df.tsne global 
tSNEmatrix$val <- df.tsne
mode$m <- 0
}

if("useheader" %in% isolate(input$clustLabels)) {
	current.clustID <- as.data.frame(Idents(scObject$val))
	setDT(current.clustID, keep.rownames = TRUE)[]
	colnames(current.clustID) <- c("clust", "ident")
	current.clustID$clust <- gsub("\\..*|_.*|-.*", "", current.clustID$clust)
	current.clustID %>% group_by(ident,clust) %>% tally() -> lables
	new.lables <- as.data.frame(lables)
	new.lables <- new.lables[order(-new.lables$n), ]
	ClustLab <- subset(new.lables, !duplicated(ident))
	new.ident <- as.vector(ClustLab$ident)
	new.clust <- as.vector(ClustLab$clust)

	#scObject$val <- RenameIdents(scObject$val, new.clust)
  Idents(scObject$val) <- plyr::mapvalues(x = Idents(scObject$val), from = new.ident, to = new.clust)
}

if(input$SwitchHeatmap == "TRUE"){
    if(is.null(scObjAllmarkers$val)){
        scObjAllmarkers$val <- FindAllMarkers(object = scObject$val, test.use = "MAST")
    }
    topN <- scObjAllmarkers$val %>% group_by(cluster) %>% top_n(input$TopDiffHeatmap, avg_logFC)
    
    DownloadPlot$val$Heatmap <- DoHeatmap(object = scObject$val, genes.use = topN$gene, rotate.key=TRUE, cex.row = 12,
                                  slim.col.label = TRUE, remove.key = FALSE)
    DownloadPlot$val$Heatmaply <- ggplotly(DownloadPlot$val$Heatmap + scale_fill_viridis(option = "viridis")) 
} else {
    if(!is.null(input$HeatmapGenes)){
        #Five visualizations of marker gene expression
        features.plot <- as.character(unlist(strsplit(input$HeatmapGenes,","))) 
        #Single cell heatmap of gene expression

        noGenes <- c()
        Genes <- c()
    
        for (i in features.plot){
          if(i %in% GetAssayData(object = scObject$val)@Dimnames[[1]])
          {
              Genes[length(Genes)+1] = i
          } else {
              noGenes[length(noGenes)+1] = i
          }
        }

        GenesAbsent$val <- noGenes

        if(length(GenesAbsent$val) == length(features.plot)) {
              plot.new()
              DownloadPlot$val$Heatmap <- NULL
              DownloadPlot$val$Heatmaply <- NULL
        } else {
              DownloadPlot$val$Heatmap <- DoHeatmap(object = scObject$val, genes.use = Genes, rotate.key=TRUE, cex.row = 12,
                                            slim.col.label = TRUE, remove.key = FALSE)
              DownloadPlot$val$Heatmaply <- ggplotly(DownloadPlot$val$Heatmap + scale_fill_viridis(option = "viridis")) 
        }
    } else {
        plot.new()
        DownloadPlot$val$Heatmap <- NULL
        DownloadPlot$val$Heatmaply <- NULL
    }
}