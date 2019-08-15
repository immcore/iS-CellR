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

########PCAplot using ggplot2
# Create data frame of clusters computed by Seurat
df.cluster <- data.frame(Cell = names(Idents(object = scObject$val)), Cluster = Idents(object = scObject$val))

# Create data frame of PCs compute by Seurat
df.pc <- data.frame(Embeddings(object = scObject$val, reduction = "pca"))
colnames(df.pc) <- gsub("_", "", colnames(df.pc))
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

df.PCAClusters %>%
  dplyr::group_by(Cluster) %>%
  dplyr::summarize(x = median(x = PC1), y = median(x = PC2)) -> labCluster

ClustID <- df.PCAClusters[tail(seq_along(df.PCAClusters),2)] # Select last 2 columns
lables <- as.data.frame(table(ClustID))
lables <- lables[order(-lables$Freq), ] # count freq of clusters with celltype 
ClustLab <- subset(lables, !duplicated(Cluster)) # remove duplicate or low counts
labCount <- as.data.frame(table(ClustID[,1])) # count freq on cluster
labCount <- labCount[order(-labCount$Freq), ]
ClustLab <- cbind(ClustLab,labCount) 
ClustInfo <- merge(labCluster, ClustLab, by = "Cluster")
ClustInfo <- ClustInfo[,c(1,2,3,4,7)]
#ClustInfo <- labCluster
names(ClustInfo) <- c("Cluster", "PC1", "PC2",  "Celltype", "Freq")
Clustd <- ClustInfo[,c(1,5)]
df.PCAClusters <- merge(df.PCAClusters, Clustd, by = "Cluster")

if("useheader" %in% isolate(input$clustLabels)) {
	df.PCAClusters$Cluster <- df.PCAClusters$Celltype
	ClustInfo$Cluster <- ClustInfo$Celltype
} 

if(input$SwitchLabelPCA == "TRUE") {
		DownloadPlot$val$PCAplot <- ggplotly(ggplot(data = df.PCAClusters, mapping = aes(PC1, PC2)) +
			geom_point(size = 1, mapping = aes(text=sprintf("Cell: %s<br>nCells: %s", Cell, Freq), color = Cluster)) +
    		#geom_point(data = ClustInfo, colour="black", size = 0, alpha = 0) +
    		geom_text(data = ClustInfo, mapping = aes(label = Cluster), size = 5, colour="black") + 
    		theme(legend.title = element_blank()),
      		tooltip = c("text"))

 	 	DownloadPlot$val$PCAplot2D <- ggplot(data = df.PCAClusters, mapping = aes(PC1, PC2)) +
			geom_point(size = 1, mapping = aes(text=sprintf("Cell: %s<br>nCells: %s", Cell, Freq), color = Cluster)) +
    		#geom_point(data = ClustInfo, colour="black", size = 0, alpha = 0) +
    		geom_text(data = ClustInfo, mapping = aes(label = Cluster), size = 5, colour="black") + 
    		theme(legend.title = element_blank())
	} else {
		DownloadPlot$val$PCAplot <- ggplotly(ggplot(data = df.PCAClusters, mapping = aes(PC1, PC2)) +
			geom_point(size = 1, mapping = aes(text=sprintf("Cell: %s<br>nCells: %s", Cell, Freq), color = Cluster)) +
			theme(legend.title = element_blank()),
    		tooltip = c("text"))

		DownloadPlot$val$PCAplot2D <- ggplot(data = df.PCAClusters, mapping = aes(PC1, PC2)) +
			geom_point(size = 1, mapping = aes(text=sprintf("Cell: %s<br>nCells: %s", Cell, Freq), color = Cluster)) +
			theme(legend.title = element_blank())
		} 
