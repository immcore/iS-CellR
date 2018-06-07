#!/usr/bin/env Rscript

if(!isS4(tSNEObj$val) || (mode$n == 0 && mode$m == 1) || "defLabels" %in% isolate(input$clustLabels) || input$changeLabels)
{
	if(!isS4(tSNEClusters$val) || (mode$n == 0 && mode$m == 1))
	{	
		Featuredata <- ProjectPCA(object = seuratObject$val, do.print = FALSE)

		Featuredata <- FindClusters(object = Featuredata, reduction.type = "pca", dims.use = 1:10, 
                            resolution = 0.6, print.output = 0, save.SNN = TRUE)

		# Make tsne global 
		#assign('tSNEClusters', Featuredata, envir=.GlobalEnv)
    tSNEClusters$val <- Featuredata
	}

	if(!isS4(tSNEObj$val) || (mode$n == 0 && mode$m == 1))
	{
		Featuredata <- RunTSNE(object = tSNEClusters$val, dims.use = 1:10, do.fast = TRUE)
		# Make tsne.obj global 
		#assign('tSNEObj', Featuredata, envir=.GlobalEnv)
    tSNEObj$val <- Featuredata
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
#assign('tSNEmatrix', df.tsne, envir=.GlobalEnv)
tSNEmatrix$val <- df.tsne
mode$m <- 0

}
#Five visualizations of marker gene expression
features.plot <- as.character(unlist(strsplit(input$FeatureGenes,","))) 

########### tSNE plot ggplot2 
# Create data frame of clusters computed by Seurat

plot_list <- list()
noGenes <- c()
df_sorted_expr <- data.frame()
df_sorted_expr_plot <- data.frame()

for (i in features.plot) { 

	if(i %in% tSNEObj$val@data@Dimnames[[1]])
	{
		df <- data.frame(Expression = tSNEObj$val@data[i,], Gene = i)
		# Add Cell column
		df$Cell <- rownames(df)

		# Merge tSNE data frame to Cluster data frame
		df <- merge(df, tSNEmatrix$val, by = "Cell")

		df_sorted_expr <- df[with(df, order(Expression)), ]

		df_sorted_expr_plot <- rbind(df_sorted_expr_plot, df_sorted_expr)
	}
	else{
		noGenes[length(noGenes)+1] = i
	}
}

#assign('GenesAbsent', noGenes, envir=.GlobalEnv)
GenesAbsent$val <- noGenes

names(df_sorted_expr_plot) <- c("Cell", "Expression", "Gene", "tSNE_1", "tSNE_2", "Cluster")
#df_sorted_expr_plot <- df_sorted_expr_plot[,c("Cell", "Expression", "Gene", "Cluster", "tSNE_1", "tSNE_2",)]

df_sorted_expr_plot$Celltype <- df_sorted_expr_plot$Cell
df_sorted_expr_plot$Celltype <- gsub("\\..*|_.*|-.*", "", df_sorted_expr_plot$Celltype)
df_sorted_expr_plot <- df_sorted_expr_plot[,c("Cell", "Expression", "Gene", "tSNE_1", "tSNE_2","Celltype","Cluster")]

df_sorted_expr_plot %>%
	dplyr::group_by(Cluster) %>%
	summarize(x = median(x = tSNE_1), y = median(x = tSNE_2)) -> labCluster

ClustID <- df_sorted_expr_plot[tail(seq_along(df_sorted_expr_plot),5)] # Select last 2 columns
ClustCell <- ClustID[,c(5,4)]
ClustGene <- ClustID[,c(5,1)]

lables <- as.data.frame(table(ClustGene))
lables <- lables[order(-lables$Freq), ] # count freq of clusters with celltype 
ClustLab <- subset(lables, !duplicated(Cluster)) # remove duplicate or low counts
ClustInfo <- merge(labCluster, ClustLab, by = "Cluster")
names(ClustInfo) <- c("Cluster", "tSNE_1", "tSNE_2",  "Gene", "Freq")
#write.table(ClustInfo, "features.txt", sep = "\t", col.names = TRUE, row.names = FALSE)
Clustd <- ClustInfo[,c(1,5)]

df_sorted <- merge(df_sorted_expr_plot, Clustd, by = "Cluster")
df_sorted_expr_plot <- df_sorted[with(df_sorted, order(Expression, Gene)), ]

ExprGenes <- as.data.frame(subset(df_sorted_expr_plot, Expression > 0), quote=FALSE)

#ExprGenes <- ExprGenes[,c(1,4,2)]

ExprGenes$allClust <- table(ExprGenes$Gene)[ExprGenes$Gene]

setDT(ExprGenes)
ExprGenes[, uniqClust := .N, by = .(Cluster, Gene)]

#write.table(ExprGenes, "features.txt", sep = "\t", col.names = TRUE, row.names = FALSE)

lables <- as.data.frame(table(ClustCell))
lables <- lables[order(-lables$Freq), ] # count freq of clusters with celltype 
ClustLab <- subset(lables, !duplicated(Cluster)) # remove duplicate or low counts
labCount <- as.data.frame(table(ClustCell[,1])) # count freq on cluster
labCount <- labCount[order(-labCount$Freq), ]
ClustLab <- cbind(ClustLab,labCount) 
ClustInfo <- merge(labCluster, ClustLab, by = "Cluster")
ClustInfo <- ClustInfo[,c(1,2,3,4,7)]
names(ClustInfo) <- c("Cluster", "tSNE_1", "tSNE_2",  "Celltype", "Freq")
#write.table(ClustInfo, "features.txt", sep = "\t", col.names = TRUE, row.names = FALSE)
#Clustd <- ClustInfo[,c(1,5)]

plot_list <- length(unique(df_sorted_expr_plot$Gene))

if(plot_list == 1){
	lab_size = 5
} else if(plot_list == 2){
	lab_size = 4
} else if(plot_list >= 3 && plot_list <= 9){
	lab_size = 3
} else if(plot_list > 9){
	lab_size = 2
}

if("useheader" %in% isolate(input$clustLabels)) {
	df_sorted_expr_plot$Cluster <- df_sorted_expr_plot$Celltype
	ClustInfo$Cluster <- ClustInfo$Celltype
  ExprGenes$Cluster <- ExprGenes$Celltype
} 

if(length(GenesAbsent) == length(features.plot)) {
	plot.new()
	} else if(input$PrintLabel_Featplot) {
	DownloadPlot$val$Featureplot <- ggplotly(ggplot(data = df_sorted_expr_plot, aes(tSNE_1, tSNE_2)) + 
      geom_point(mapping = aes(text=sprintf("Cell: %s<br>Cluster: %s [nCells: %s]", Cell, Cluster, Freq), color = Expression)) +
      geom_point(data=ExprGenes, mapping=aes(text=sprintf("Cell: %s<br>Cluster: %s [nCells: %s/%s]", Cell, Cluster, uniqClust, Freq), color = Expression)) +
      scale_color_gradientn(limits = c(min(df_sorted_expr_plot$Expression),max(df_sorted_expr_plot$Expression)), colors= c("lightgrey", "blue", "red"), name = i, guide = guide_legend(title.hjust = 1)) + 
      ggtitle("") + theme(plot.title = element_text(color="black", size=14, face="bold")) + theme(legend.position="right") + theme(legend.title = element_text()) +
      facet_wrap(~Gene, scales = "free") + 
      geom_text(data = ClustInfo, mapping = aes(label = Cluster), size = lab_size, colour="black") +
      theme(strip.background = element_blank(), strip.placement = "outside") +
      theme(strip.text = element_text(size = rel(1))) +
      xlab("tSNE1") + ylab("tSNE2") +
      theme(axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 2, 0, 0), "mm"), angle = 90)) +
      theme(plot.margin = unit(c(1,1,1,1), "cm"), panel.margin = unit(0, "cm")),
      #theme(panel.spacing=unit(4, "cm")) +
      tooltip = c("text"))
   	DownloadPlot$val$Featureplot2D <- ggplot(data = df_sorted_expr_plot, aes(tSNE_1, tSNE_2)) + 
      geom_point(mapping = aes(text=sprintf("Cell: %s<br>Cluster: %s [nCells: %s]", Cell, Cluster, Freq), color = Expression)) +
      geom_point(data=ExprGenes, mapping=aes(text=sprintf("Cell: %s<br>Cluster: %s [nCells: %s/%s]", Cell, Cluster, uniqClust, Freq), color = Expression)) +
      scale_color_gradientn(limits = c(min(df_sorted_expr_plot$Expression),max(df_sorted_expr_plot$Expression)), colors= c("lightgrey", "blue", "red"), name = i, guide = guide_legend(title.hjust = 1)) + 
      ggtitle("") + theme(plot.title = element_text(color="black", size=14, face="bold")) + theme(legend.position="right") + theme(legend.title = element_text()) +
      facet_wrap(~Gene, scales = "free") + 
      geom_text(data = ClustInfo, mapping = aes(label = Cluster), size = lab_size, colour="black") +
      theme(strip.background = element_blank(), strip.placement = "outside") +
      theme(strip.text = element_text(size = rel(1))) +
      xlab("tSNE1") + ylab("tSNE2") +
      theme(axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 2, 0, 0), "mm"), angle = 90)) +
      theme(plot.margin = unit(c(1,1,1,1), "cm"), panel.margin = unit(0, "cm"))
      #theme(panel.spacing=unit(4, "cm")) +
      } else {
	DownloadPlot$val$Featureplot <- ggplotly(ggplot(data = df_sorted_expr_plot, aes(tSNE_1, tSNE_2)) + 
      geom_point(mapping = aes(text=sprintf("Cell: %s<br>Cluster: %s [nCells: %s]", Cell, Cluster, Freq), color = Expression)) +
      scale_color_gradientn(limits = c(min(df_sorted_expr_plot$Expression),max(df_sorted_expr_plot$Expression)), colors= c("lightgrey", "blue", "red"), name = i, guide = guide_legend(title.hjust = 1)) + 
      ggtitle("") + theme(plot.title = element_text(color="black", size=14, face="bold")) + theme(legend.position="right") + theme(legend.title = element_text()) +
      facet_wrap(~Gene, scales = "free") + 
      geom_point(data=ExprGenes, mapping=aes(text=sprintf("Cell: %s<br>Cluster: %s [nCells: %s/%s]", Cell, Cluster, uniqClust, Freq),color = Expression)) +
      theme(strip.background = element_blank(), strip.placement = "outside") +
      theme(strip.text = element_text(size = rel(1))) +
      xlab("tSNE1") + ylab("tSNE2") +
      theme(axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")),
      axis.title.y = element_text(margin = unit(c(0, 2, 0, 0), "mm"), angle = 90)) +
      theme(plot.margin = unit(c(1,1,1,1), "cm"), panel.margin = unit(0, "cm")),
      #theme(panel.spacing=unit(4, "cm")) +
      tooltip = c("text"))
   	DownloadPlot$val$Featureplot2D <- ggplot(data = df_sorted_expr_plot, aes(tSNE_1, tSNE_2)) + 
      geom_point(mapping = aes(text=sprintf("Cell: %s<br>Cluster: %s [nCells: %s]", Cell, Cluster, Freq), color = Expression)) +
      geom_point(data=ExprGenes, mapping=aes(text=sprintf("Cell: %s<br>Cluster: %s [nCells: %s/%s]", Cell, Cluster, uniqClust, Freq), color = Expression)) +
      scale_color_gradientn(limits = c(min(df_sorted_expr_plot$Expression),max(df_sorted_expr_plot$Expression)), colors= c("lightgrey", "blue", "red"), name = i, guide = guide_legend(title.hjust = 1)) + 
      ggtitle("") + theme(plot.title = element_text(color="black", size=14, face="bold")) + theme(legend.position="right") + theme(legend.title = element_text()) +
      facet_wrap(~Gene, scales = "free") + 
      theme(strip.background = element_blank(), strip.placement = "outside") +
      theme(strip.text = element_text(size = rel(1))) +
      xlab("tSNE1") + ylab("tSNE2") +
      theme(axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")),
      axis.title.y = element_text(margin = unit(c(0, 2, 0, 0), "mm"), angle = 90)) +
      theme(plot.margin = unit(c(1,1,1,1), "cm"), panel.margin = unit(0, "cm"))
      #theme(panel.spacing=unit(4, "cm")) +
   	} 



	

