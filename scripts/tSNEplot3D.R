#!/usr/bin/env Rscript

if(!isS4(tSNEClusters$val) || "defLabels"  %in% isolate(input$clustLabels) || input$changeLabels)
{
    if(!isS4(tSNEClusters$val)){
        tSNEdata <- ProjectPCA(object = seuratObject$val, do.print = FALSE)
        tSNEdata <- FindClusters(object = tSNEdata, reduction.type = "pca", dims.use = 1:10, 
                            resolution = 0.6, print.output = 0, save.SNN = TRUE)
        # Make df.tsne global 
        tSNEClusters$val <- tSNEdata
    }#mode$n <- 1

if(!isS4(tSNE3D$val))
    {
        tSNE3D$val <- RunTSNE(object = tSNEClusters$val, dims.use = 1:10, dim.embed= 3, do.fast = TRUE)
    }

    ################## for Custom labels ################
  if("customLabels" %in% isolate(input$clustLabels)) {
    cluster.ids <- as.character(unlist(ClusterLabInfo$val[,1]))#, decreasing = FALSE)#c("CD4", "Bcells", "CD8cells",
    new.cluster.ids <- as.character(unlist(ClusterLabInfo$val[,2]))#c("CD4", "Bcells", "CD8cells", 
    tSNE3D$val@ident <- plyr::mapvalues(x = tSNE3D$val@ident, from = cluster.ids, to = new.cluster.ids)
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
df.tsne <- data.frame(tSNE3D$val@dr$tsne@cell.embeddings)
# Add Cell column
df.tsne$Cell = rownames(df.tsne)
# Merge tSNE data frame to Cluster data frame
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

DownloadPlot$val$tSNE3D <- plot_ly(df.tsne, x = ~tSNE_1, y = ~tSNE_2, z = ~tSNE_3, type = "scatter3d", mode = "markers", color = ~Cluster, 
            text = ~paste("Cell:", Cell, "<br>Cluster:", Cluster)) 
