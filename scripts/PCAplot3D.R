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

########PCAplot using ggplot2
# Create data frame of clusters computed by Seurat
df.cluster <- data.frame(Cell = names(scObject$val@ident), Cluster = scObject$val@ident)

# Create data frame of PCs compute by Seurat
df.pc <- data.frame(scObject$val@dr$pca@cell.embeddings)
# Add Cell column
df.pc$Cell = rownames(df.pc)
# Merge PC data frame to Cluster data frame
df.pc <- merge(df.pc, df.cluster, by = "Cell")
# Plot PC1 and PC2 using ggplot2,
# colors of points correspond to Cluster IDs
# Make df.tsne global 

#assign('PCAClusters', df.pc, envir=.GlobalEnv)
PCAClusters$val <- df.pc
mode$l <- 0
}

df.PCAClusters <- PCAClusters$val

df.PCAClusters$Celltype <- df.PCAClusters$Cell
df.PCAClusters$Celltype <- gsub("\\..*|_.*|-.*", "", df.PCAClusters$Celltype)

if("useheader" %in% isolate(input$clustLabels)) {
    df.PCAClusters$Cluster <- df.PCAClusters$Celltype
} 

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

DownloadPlot$val$PCA3D <- plot_ly(df.PCAClusters, x = ~PC1, y = ~PC2, z = ~PC3, type = "scatter3d", mode = "markers", color = ~Cluster, 
    		text = ~paste("Cell: ", Cell, "\nCluster: ", Cluster)) 
    
