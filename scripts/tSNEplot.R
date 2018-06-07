#!/usr/bin/env Rscript

if(!isS4(tSNEObj$val) || (mode$n == 0 && mode$m == 1) || "defLabels" %in% isolate(input$clustLabels) || input$changeLabels)
{
	if(!isS4(tSNEClusters$val) || (mode$n == 0 && mode$m == 1))
	{
		tSNEdata <- ProjectPCA(object = seuratObject$val, do.print = FALSE)
 
		tSNEdata <- FindClusters(object = tSNEdata, reduction.type = "pca", dims.use = 1:10, 
                            resolution = 0.6, print.output = 0, save.SNN = TRUE)

		# Make df.tsne global 
		tSNEClusters$val <- tSNEdata
	}

	if(!isS4(tSNEObj$val) || (mode$n == 0 && mode$m == 1))
	{
		tSNEdata <- RunTSNE(object = tSNEClusters$val, dims.use = 1:10, do.fast = TRUE)
  
		# Make tsne.obj global 
		tSNEObj$val <- tSNEdata
	}

	################## for Custom labels ################
  if("customLabels" %in% isolate(input$clustLabels)) {
    cluster.ids <- as.character(unlist(ClusterLabInfo$val[,1]))#, decreasing = FALSE)#c("CD4", "Bcells", "CD8cells",
    new.cluster.ids <- as.character(unlist(ClusterLabInfo$val[,2]))#c("CD4", "Bcells", "CD8cells", 
    tSNEObj$val@ident <- plyr::mapvalues(x = tSNEObj$val@ident, from = cluster.ids, to = new.cluster.ids)
    CInfo <- cbind(cluster.ids,new.cluster.ids)
    ClusterLabInfo$val <- CInfo
  }

  if("defLabels" %in% isolate(input$clustLabels)) {
    if(!is.null(dfcluster.ids$val)){
      new.cluster.ids <- as.character(unlist(ClusterLabInfo$val[,2]))
      current.ids <- as.character(unlist(ClusterLabInfo$val[,1]))#, decreasing = FALSE)#c("CD4", "Bcells", "CD8cells",
      tSNEObj$val@ident <- plyr::mapvalues(x = tSNEObj$val@ident, from = new.cluster.ids, to = current.ids)
    } else {
      new.cluster.ids = ""
      current.ids <- sort(as.character(unique(tSNEObj$val@ident)), decreasing = FALSE)
    }
    cluster.ids <- current.ids
    CInfo <- cbind(cluster.ids,new.cluster.ids)
    ClusterLabInfo$val <- CInfo 
    dfcluster.ids$val <- cluster.ids      
  }

########### tSNE plot ggplot2 
# Create data frame of clusters computed by Seurat
df.cluster <- data.frame(Cell = names(tSNEObj$val@ident), Cluster = tSNEObj$val@ident)
    
# Create data frame of tSNE compute by Seurat
df.tsne <- data.frame(tSNEObj$val@dr$tsne@cell.embeddings)
# Add Cell column
df.tsne$Cell = rownames(df.tsne)
# Merge tSNE data frame to Cluster data frame
df.tsne <- merge(df.tsne, df.cluster, by = "Cell")

# Make df.tsne global 
tSNEmatrix$val <- df.tsne
mode$m <- 0
}

df.tSNEmatrix <- tSNEmatrix$val

df.tSNEmatrix$Celltype <- df.tSNEmatrix$Cell
df.tSNEmatrix$Celltype <- gsub("\\..*|_.*|-.*", "", df.tSNEmatrix$Celltype)

df.tSNEmatrix %>%
  dplyr::group_by(Cluster) %>%
  summarize(x = median(x = tSNE_1), y = median(x = tSNE_2)) -> labCluster

ClustID <- df.tSNEmatrix[tail(seq_along(df.tSNEmatrix),2)] # Select last 2 columns
#write.table(ClustID, "features.txt", sep = "\t", col.names = TRUE, row.names = FALSE)
  
lables <- as.data.frame(table(ClustID))
lables <- lables[order(-lables$Freq), ] # count freq of clusters with celltype 
ClustLab <- subset(lables, !duplicated(Cluster)) # remove duplicate or low counts
labCount <- as.data.frame(table(ClustID[,1])) # count freq on cluster
labCount <- labCount[order(-labCount$Freq), ]
ClustLab <- cbind(ClustLab,labCount) 
ClustInfo <- merge(labCluster, ClustLab, by = "Cluster")
ClustInfo <- ClustInfo[,c(1,2,3,4,7)]

#ClustInfo <- labCluster
names(ClustInfo) <- c("Cluster", "tSNE_1", "tSNE_2",  "Celltype", "Freq")
Clustd <- ClustInfo[,c(1,5)]
df.tSNEmatrix <- merge(df.tSNEmatrix, Clustd, by = "Cluster")

if("useheader" %in% isolate(input$clustLabels)) {
	df.tSNEmatrix$Cluster <- df.tSNEmatrix$Celltype
	ClustInfo$Cluster <- ClustInfo$Celltype
} 

if(input$PrintLabeltSNE) {
		DownloadPlot$val$tSNEplot <- ggplotly(ggplot(data = df.tSNEmatrix, aes(tSNE_1, tSNE_2)) + 
			geom_point(mapping = aes(text=sprintf("Cell: %s<br>nCells: %s", Cell, Freq), color = Cluster)) +
            #geom_point(data = labCluster, mapping = aes(x = labCluster$x, y = labCluster$y), colour="black", size = 0, alpha = 0) +
            geom_text(data = ClustInfo, mapping = aes(label = Cluster), size = 5, colour="black") +
            theme(legend.title = element_blank()),
            tooltip = c("text"))

 	 	DownloadPlot$val$tSNEplot2D <- ggplot(data = df.tSNEmatrix, aes(tSNE_1, tSNE_2)) + 
			geom_point(mapping = aes(text=sprintf("Cell: %s<br>nCells: %s", Cell, Freq), color = Cluster)) +
            #geom_point(data = labCluster, mapping = aes(x = labCluster$x, y = labCluster$y), colour="black", size = 0, alpha = 0) +
            geom_text(data = ClustInfo, mapping = aes(label = Cluster), size = 5, colour="black") +
            theme(legend.title = element_blank())

	} else {
		DownloadPlot$val$tSNEplot <- ggplotly(ggplot(data = df.tSNEmatrix, aes(tSNE_1, tSNE_2)) + 
			geom_point(mapping = aes(text=sprintf("Cell: %s<br>nCells: %s", Cell, Freq), color = Cluster)) +
			theme(legend.title = element_blank()),
			tooltip = c("text"))

		DownloadPlot$val$tSNEplot2D <- ggplot(data = df.tSNEmatrix, aes(tSNE_1, tSNE_2)) + 
			geom_point(mapping = aes(text=sprintf("Cell: %s<br>nCells: %s", Cell, Freq), color = Cluster)) +
			theme(legend.title = element_blank())
		}
 