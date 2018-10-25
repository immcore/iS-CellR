#!/usr/bin/env Rscript

if(isS4(scObject$val) || "defLabels"  %in% isolate(input$clustLabels) || "useheader" %in% isolate(input$clustLabels) || input$changeLabels)
{
  if(!isS4(tSNE3D$val))
    {
        tSNE3D$val <- RunTSNE(object = scObject$val, reduction.use = "pca", dims.use = 1:75, dim.embed = 3, tsne.method = "FIt-SNE", nthreads = 4, reduction.name = "FItSNE", reduction.key = "FItSNE_", fast_tsne_path = "/bin/fast_tsne", max_iter = 2000)
        tSNE3D$val <- RunUMAP(object = tSNE3D$val, reduction.use = "pca", max.dim = 3, dims.use = 1:75, min_dist = 0.75)
    }

    ################## for Custom labels ################
  if("customLabels" %in% isolate(input$clustLabels)) {
    cluster.ids <- as.character(unlist(ClusterLabInfo$val[,1]))#, decreasing = FALSE)#c("CD4", "Bcells", "CD8cells",
    new.cluster.ids <- as.character(unlist(ClusterLabInfo$val[,2]))#c("CD4", "Bcells", "CD8cells", 
    if(input$changeLabels){
      tSNE3D$val@ident <- plyr::mapvalues(x = tSNE3D$val@ident, from = cluster.ids, to = new.cluster.ids)
    }
    CInfo <- cbind(cluster.ids,new.cluster.ids)
    ClusterLabInfo$val <- CInfo
  }

  if("defLabels" %in% isolate(input$clustLabels)) {
    if(!is.null(dfcluster.ids$val)){
      new.cluster.ids <- as.character(unlist(ClusterLabInfo$val[,2]))
      current.ids <- as.character(unlist(ClusterLabInfo$val[,1]))#, decreasing = FALSE)#c("CD4", "Bcells", "CD8cells",
      tSNE3D$val@ident <- plyr::mapvalues(x = tSNE3D$val@ident, from = new.cluster.ids, to = current.ids)
    } else {
      new.cluster.ids = ""
      current.ids <- sort(as.character(unique(tSNE3D$val@ident)), decreasing = FALSE)
    }
    cluster.ids <- current.ids
    CInfo <- cbind(cluster.ids,new.cluster.ids)
    ClusterLabInfo$val <- CInfo 
    dfcluster.ids$val <- cluster.ids      
  }

}
 
########### tSNE plot ggplot2
# Create data frame of clusters computed by Seurat
df.cluster <- data.frame(Cell = names(tSNE3D$val@ident), Cluster = tSNE3D$val@ident)
    
# Create data frame of tSNE compute by Seurat
df.umap <- data.frame(tSNE3D$val@dr$umap@cell.embeddings)
# Add Cell column
df.umap$Cell = rownames(df.umap)
# Create data frame of tSNE compute by Seurat
df.FItsne <- data.frame(tSNE3D$val@dr$FItSNE@cell.embeddings)
# Add Cell column
df.FItsne$Cell = rownames(df.FItsne)

# Merge tSNE data frame to Cluster data frame
df.tsne <- merge(df.umap, df.FItsne, by = "Cell")
df.tsne <- merge(df.tsne, df.cluster, by = "Cell")

df.tsne$Celltype <- df.tsne$Cell
#df.tsne$Celltype <- gsub("CY|cy|Cy|cY", "Mel", df.tsne$Celltype)
df.tsne$Celltype <- gsub("\\..*|_.*|-.*", "", df.tsne$Celltype)

xaxis <- list(
    #title = names(x = data.plot)[1],
    showgrid = TRUE,
    zeroline = FALSE,
    showticklabels = TRUE,
    showline = TRUE
  )
yaxis <- list(
    #title = names(x = data.plot)[2],
    showgrid = TRUE,
    zeroline = FALSE,
    showticklabels = TRUE,
    showline = TRUE
  )
zaxis <- list(
    #title = names(x = data.plot)[2],
    showgrid = TRUE,
    zeroline = FALSE,
    showticklabels = TRUE,
    showline = TRUE
  )

if("useheader" %in% isolate(input$clustLabels)) {
    df.tsne$Cluster <- df.tsne$Celltype
} 

if("UMAP" %in% isolate(dimPkg$val)) {
  DownloadPlot$val$tSNE3D <- plot_ly(df.tsne, x = ~UMAP1, y = ~UMAP2, z = ~UMAP3, type = "scatter3d", mode = "markers", color = ~Cluster, 
            text = ~paste("Cell:", Cell, "<br>Cluster:", Cluster)) 
} else {
  DownloadPlot$val$tSNE3D <- plot_ly(df.tsne, x = ~FItSNE_1, y = ~FItSNE_2, z = ~FItSNE_3, type = "scatter3d", mode = "markers", color = ~Cluster, 
            text = ~paste("Cell:", Cell, "<br>Cluster:", Cluster)) 
}
