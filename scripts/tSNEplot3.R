#!/usr/bin/env Rscript

# Plot tSNE1 and tSNE2 using ggplot2,
if(!isS4(tSNEObj$val) || (mode$n == 0 && mode$m == 1) || "defLabels"  %in% isolate(input$clustLabels) || input$changeLabels)
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
lables <- as.data.frame(table(ClustID))
lables <- lables[order(-lables$Freq), ]
ClustLab <- subset(lables, !duplicated(Cluster))
ClustInfo <- merge(labCluster, ClustLab, by = "Cluster")
names(ClustInfo) <- c("Cluster", "tSNE_1", "tSNE_2",  "Celltype", "Freq")

if("useheader" %in% isolate(input$clustLabels)) {
	df.tSNEmatrix$Cluster <- df.tSNEmatrix$Celltype
	ClustInfo$Cluster <- ClustInfo$Celltype
}

s1 <- input$tSNEplot3_rows_current  # rows on the current page
s2 <- input$tSNEplot3_rows_all      # rows on all pages (after being filtered)
s3 <- input$tSNEplot3_rows_selected  # rows on the current page

p1 <- ggplot(data = df.tSNEmatrix, mapping=aes(tSNE_1, tSNE_2)) + geom_point(aes(color = Cluster)) + theme(legend.title = element_blank()) 

print(p1)

if (length(s2) > 0 && length(s2) < nrow(df.tSNEmatrix)) {
	# show red circles when performing searching
	#p1 + geom_point(data=tSNEmatrix$val[s2, , drop = FALSE], mapping=aes(tSNE_1, tSNE_2, color=Cluster), size=7)
	p1 + geom_point(data=df.tSNEmatrix[s2, , drop = FALSE], mapping=aes(tSNE_1, tSNE_2), colour="black", size=2)
	#if (length(s1)){p1 + geom_point(data=tSNEmatrix$val[s3, , drop = FALSE], mapping=aes(tSNE_1, tSNE_2), colour="red", shape=21, stroke = 2, size=5)}
} else if (length(s1)) {
	# solid dots (pch = 19) for current page
	p1 + geom_point(data=df.tSNEmatrix[s1, , drop = FALSE], mapping=aes(tSNE_1, tSNE_2), colour="black", size=2)
	#if (length(s1)){p1 + geom_point(data=tSNEmatrix$val[s3, , drop = FALSE], mapping=aes(tSNE_1, tSNE_2), colour="red", shape=21, stroke = 2, size=5)}
} else {
	print(p1)
}

