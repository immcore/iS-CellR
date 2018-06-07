#!/usr/bin/env Rscript

if(!isS4(tSNEObj$val) || (mode$n == 0 && mode$m == 1) || "defLabels" %in% isolate(input$clustLabels) || input$changeLabels)
{
  if(!isS4(tSNEClusters$val) || (mode$n == 0 && mode$m == 1))
  {
    CoExprdata <- ProjectPCA(object = seuratObject$val, do.print = FALSE)
     
    CoExprdata <- FindClusters(object = CoExprdata, reduction.type = "pca", dims.use = 1:10, 
                            resolution = 0.6, print.output = 0, save.SNN = TRUE)

    # Make df.tsne global 
    tSNEClusters$val <- CoExprdata
  }
  
  if(!isS4(tSNEObj$val) || (mode$n == 0 && mode$m == 1))
  {
    CoExprdata <- RunTSNE(object = tSNEClusters$val, dims.use = 1:10, do.fast = TRUE)
    # Make df.tsne global 
    tSNEObj$val <- CoExprdata
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

CoExprPLOT <- function(){
  CoExpr <- read.delim("CoExpr.txt", sep="\t", header=F)
names(CoExpr) <- c("Cell", "Expression", "Gene", "tSNE_1", "tSNE_2", "Cluster")
# Plot gene expression using ggplot2,
# Colors of points correspond to gene expression levels

CoExpr$Celltype <- CoExpr$Cell
CoExpr$Celltype <- gsub("\\..*|_.*|-.*", "", CoExpr$Celltype)

CoExprValue$val$Gene1 <- features.plot[1]
CoExprValue$val$Gene2 <- features.plot[2]

CoExprValue$val$Gene1Min <- min(subset(CoExpr, Gene == features.plot[1])$Expression)
CoExprValue$val$Gene1Max <- max(subset(CoExpr, Gene == features.plot[1])$Expression)
CoExprValue$val$Gene1Mean <- mean(subset(CoExpr, Gene == features.plot[1])$Expression)

CoExprValue$val$Gene2Min <- min(subset(CoExpr, Gene == features.plot[2])$Expression)
CoExprValue$val$Gene2Max <- max(subset(CoExpr, Gene == features.plot[2])$Expression)
CoExprValue$val$Gene2Mean <- mean(subset(CoExpr, Gene == features.plot[2])$Expression)

Gene1 <- as.data.frame(subset(CoExpr, Gene == features.plot[1] & Expression >= input$ExprCutoff))
Gene2 <- as.data.frame(subset(CoExpr, Gene == features.plot[2] & Expression >= input$ExprCutoff))

Gene12 <- merge(Gene1, Gene2, by = "Cell")
Gene12 <- Gene12[,1:9]
names(Gene12) <- c("Cell", "Expression.x", "Gene.x", "tSNE_1", "tSNE_2", "Cluster", "Celltype", "Expression.y", "Gene.y")

diff.Gene1 <- (!Gene1[,1] %in% Gene12[,1])
Gene1only <- Gene1[diff.Gene1,]
Gene1 <- Gene1only

diff.Gene2 <- (!Gene2[,1] %in% Gene12[,1])
Gene2only <- Gene2[diff.Gene2,]
Gene2 <- Gene2only

G1 <- paste0(features.plot[1]," High")
G2 <- paste0(features.plot[2]," High")
G12 <- paste0("Both High")

Gene12_cellCount <- function() {
  Gene12 %>%
  dplyr::group_by(Cluster) %>%
  summarize(x = median(x = tSNE_1), y = median(x = tSNE_2)) -> labCluster
  
  #ClustID <- CoExpr[tail(seq_along(CoExpr),2)] # Select last 2 columns
  ClustID <- Gene12[tail(seq_along(Gene12),4)] # Select last 2 columns
  ClustID <- ClustID[,c(1,2)]
  lables <- as.data.frame(table(ClustID))
  lables <- lables[order(-lables$Freq), ] # count freq of clusters with celltype 
  ClustLab <- subset(lables, !duplicated(Cluster)) # remove duplicate or low counts
  labCount <- as.data.frame(table(ClustID[,1])) # count freq on cluster
  labCount <- labCount[order(-labCount$Freq), ]
  ClustLab <- cbind(ClustLab,labCount) 
  ClustInfo <- merge(labCluster, ClustLab, by = "Cluster")
  ClustInfo <- ClustInfo[,c(1,2,3,4,7)]
  names(ClustInfo) <- c("Cluster", "tSNE_1", "tSNE_2",  "Celltype", "Freq")
  Clustd <- ClustInfo[,c(1,5)]
  df.Gene12 <- merge(Gene12, Clustd, by = "Cluster")
  #write.table(df.Gene12, "features12.txt", sep = "\t", col.names = TRUE, row.names = FALSE)
  return(df.Gene12)
}

Gene1_cellCount <- function() {
  Gene1 %>%
  dplyr::group_by(Cluster) %>%
  summarize(x = median(x = tSNE_1), y = median(x = tSNE_2)) -> labCluster

  ClustID <- Gene1[tail(seq_along(Gene1),2)] # Select last 2 columns
  lables <- as.data.frame(table(ClustID))
  lables <- lables[order(-lables$Freq), ] # count freq of clusters with celltype 
  
  ClustLab <- subset(lables, !duplicated(Cluster)) # remove duplicate or low counts
  labCount <- as.data.frame(table(ClustID[,1])) # count freq on cluster
  labCount <- labCount[order(-labCount$Freq), ]
  ClustLab <- cbind(ClustLab,labCount) 
  ClustInfo <- merge(labCluster, ClustLab, by = "Cluster")
  ClustInfo <- ClustInfo[,c(1,2,3,4,7)]
  names(ClustInfo) <- c("Cluster", "tSNE_1", "tSNE_2",  "Celltype", "Freq")
  Clustd <- ClustInfo[,c(1,5)]
  df.Gene1 <- merge(Gene1, Clustd, by = "Cluster")
  #write.table(df.Gene1, "features1.txt", sep = "\t", col.names = TRUE, row.names = FALSE)
  return(df.Gene1)
  #write.table(df.Gene1, "features1.txt", sep = "\t", col.names = TRUE, row.names = FALSE)
}

Gene2_cellCount <- function(){
  Gene2 %>%
  dplyr::group_by(Cluster) %>%
  summarize(x = median(x = tSNE_1), y = median(x = tSNE_2)) -> labCluster

  ClustID <- Gene2[tail(seq_along(Gene2),2)] # Select last 2 columns
  lables <- as.data.frame(table(ClustID))
  lables <- lables[order(-lables$Freq), ] # count freq of clusters with celltype 
  ClustLab <- subset(lables, !duplicated(Cluster)) # remove duplicate or low counts
  labCount <- as.data.frame(table(ClustID[,1])) # count freq on cluster
  labCount <- labCount[order(-labCount$Freq), ]
  ClustLab <- cbind(ClustLab,labCount) 
  ClustInfo <- merge(labCluster, ClustLab, by = "Cluster")
  ClustInfo <- ClustInfo[,c(1,2,3,4,7)]
  names(ClustInfo) <- c("Cluster", "tSNE_1", "tSNE_2",  "Celltype", "Freq")
  Clustd <- ClustInfo[,c(1,5)]
  df.Gene2 <- merge(Gene2, Clustd, by = "Cluster")
  #write.table(df.Gene2, "features2.txt", sep = "\t", col.names = TRUE, row.names = FALSE)
  return(df.Gene2)
}

CoExpr %>%
  dplyr::group_by(Cluster) %>%
  summarize(x = median(x = tSNE_1), y = median(x = tSNE_2)) -> labCluster

  #ClustID <- CoExpr[tail(seq_along(CoExpr),2)] # Select last 2 columns
  ClustID <- CoExpr[tail(seq_along(CoExpr),2)] # Select last 2 columns
  lables <- as.data.frame(table(ClustID))
  lables <- lables[order(-lables$Freq), ] # count freq of clusters with celltype 
  ClustLab <- subset(lables, !duplicated(Cluster)) # remove duplicate or low counts 
  labCount <- as.data.frame(table(ClustID[,1])) # count freq on cluster
  labCount <- labCount[order(-labCount$Freq), ]
  ClustLab <- cbind(ClustLab,labCount) 
  ClustInfo <- merge(labCluster, ClustLab, by = "Cluster")
  ClustInfo <- ClustInfo[,c(1,2,3,4,7)]
  names(ClustInfo) <- c("Cluster", "tSNE_1", "tSNE_2",  "Celltype", "Freq")
  Clustd <- ClustInfo[,c(1,5)]
  Clustd$Freq <- Clustd$Freq/2
  CoExpr <- merge(CoExpr, Clustd, by = "Cluster")
  #write.table(Clustd, "features.txt", sep = "\t", col.names = TRUE, row.names = FALSE)

if(nrow(Gene12) > 0){
  df.Gene12 = Gene12_cellCount()
} else {
  df.Gene12 <- NULL
}

if(nrow(Gene1) > 0){
  df.Gene1  = Gene1_cellCount()
} else {
  df.Gene1 <- NULL
}

if(nrow(Gene2) > 0){
  df.Gene2  = Gene2_cellCount()
} else {
  df.Gene2 <- NULL
}

if("useheader" %in% isolate(input$clustLabels)) {
  CoExpr$Cluster <- CoExpr$Celltype
  df.Gene1$Cluster <- df.Gene1$Celltype 
  df.Gene2$Cluster <- df.Gene2$Celltype
  df.Gene12$Cluster <- df.Gene12$Celltype
  ClustInfo$Cluster <- ClustInfo$Celltype
} 

  if(input$PrintLabel) {
        DownloadPlot$val$CoExprplot <- ggplot(data = CoExpr, aes(tSNE_1, tSNE_2)) + 
           geom_point(size=1, colour = "grey", aes(text=sprintf("Cell: %s<br>Cluster: %s [nCells: %s]<br>Gene: %s<br>Expr: %s", Cell, Cluster, Freq, Gene, Expression))) +
           geom_text(data=ClustInfo, mapping = aes(label = Cluster), size = 5, colour="black") +
           theme(legend.title = element_blank()) +
           ggtitle(paste0(features.plot[1]," x ",features.plot[2]," ","[Expr threshold: ",input$ExprCutoff,"]")) + 
           theme(plot.title = element_text(color="black", size=14, face="bold"))

        DownloadPlot$val$CoExprplot1 <- ggplot(data = CoExpr, aes(tSNE_1, tSNE_2)) + 
           geom_point(size=1, colour = "grey", aes(text=sprintf("Cell: %s<br>Cluster: %s [nCells: %s]<br>Gene: %s<br>Expr: %s", Cell, Cluster, Freq, Gene, Expression))) +
           geom_text(data=ClustInfo, mapping = aes(label = Cluster), size = 5, colour="black") +
           theme(legend.title = element_blank()) +
           ggtitle(paste0(features.plot[1]," x ",features.plot[2]," ","[Expr threshold: ",input$ExprCutoff,"]")) + 
           theme(plot.title = element_text(color="black", size=14, face="bold"))
    } else {
        DownloadPlot$val$CoExprplot <- ggplot(data = CoExpr, aes(tSNE_1, tSNE_2)) + 
           geom_point(size=1, colour = "grey", aes(text=sprintf("Cell: %s<br>Cluster: %s [nCells: %s]<br>Gene: %s<br>Expr: %s", Cell, Cluster, Freq, Gene, Expression))) +
           theme(legend.title = element_blank()) +
           ggtitle(paste0(features.plot[1]," x ",features.plot[2]," ","[Expr threshold: ",input$ExprCutoff,"]")) + 
           theme(plot.title = element_text(color="black", size=14, face="bold"))

        DownloadPlot$val$CoExprplot1 <- ggplot(data = CoExpr, aes(tSNE_1, tSNE_2)) + 
           geom_point(size=1, colour = "grey", aes(text=sprintf("Cell: %s<br>Cluster: %s [nCells: %s]<br>Gene: %s<br>Expr: %s", Cell, Cluster, Freq, Gene, Expression))) +
           theme(legend.title = element_blank()) +
           ggtitle(paste0(features.plot[1]," x ",features.plot[2]," ","[Expr threshold: ",input$ExprCutoff,"]")) + 
           theme(plot.title = element_text(color="black", size=14, face="bold"))
    }

  if(nrow(Gene1) > 0){
    DownloadPlot$val$CoExprplot <- DownloadPlot$val$CoExprplot + 
      geom_point(data=df.Gene1, mapping=aes(text=sprintf("Cell: %s<br>Cluster: %s [nCells: %s]<br>Gene: %s<br>Expr: %s", Cell, Cluster, Freq, Gene, Expression), colour=G1),  size=1)
    DownloadPlot$val$CoExprplot1 <- DownloadPlot$val$CoExprplot1 + 
      geom_point(data=df.Gene1, mapping=aes(colour=G1), size=1)

  } else {
    CoExprValue$val$Gene1Ab <- unlist(features.plot[1])
  } 
  if(nrow(Gene2) > 0){
    DownloadPlot$val$CoExprplot <- DownloadPlot$val$CoExprplot + 
      geom_point(data=df.Gene2, mapping=aes(text=sprintf("Cell: %s<br>Cluster: %s [nCells: %s]<br>Gene: %s<br>Expr: %s", Cell, Cluster, Freq, Gene, Expression), colour=G2), size=1)
    DownloadPlot$val$CoExprplot1 <- DownloadPlot$val$CoExprplot1 + 
      geom_point(data=df.Gene2, mapping=aes(colour=G2), size=1)  
  } else {
    CoExprValue$val$Gene2Ab <- unlist(features.plot[2])
  }  
  if(nrow(Gene12) > 0){
    DownloadPlot$val$CoExprplot <- DownloadPlot$val$CoExprplot + 
      geom_point(data=df.Gene12, mapping=aes(text=sprintf("Cell: %s<br>Cluster: %s [nCells: %s]<br>Gene: %s<br>Expr: %s<br>Gene: %s<br>Expr: %s", Cell, Cluster, Freq, Gene.x, Expression.x, Gene.y, Expression.y), colour=G12), size=1)
    DownloadPlot$val$CoExprplot1 <- DownloadPlot$val$CoExprplot1 + 
      geom_point(data=df.Gene12, mapping=aes(colour=G12), size=1)
  } else {
    CoExprValue$val$Genes <- unlist(features.plot)
  }

  if(nrow(Gene1) > 0 || nrow(Gene2) > 0 || nrow(Gene12) > 0){
    DownloadPlot$val$CoExprplot <- ggplotly(DownloadPlot$val$CoExprplot, tooltip = c("text"))
    DownloadPlot$val$CoExprplot1 <- DownloadPlot$val$CoExprplot1
  } 
  
  if (file.exists("CoExpr.txt")) file.remove("CoExpr.txt")
}

# CoExpr gene of interest
features.plot <- as.character(unlist(strsplit(input$CoExprGenes," "))) 

DownloadPlot$val$CoExprplot <- NULL
CoExprValue$val$Gene1Ab <- NULL
CoExprValue$val$Gene2Ab <- NULL

noGenes <- c()
p = 0

if(features.plot[1] %in% tSNEObj$val@data@Dimnames[[1]] && features.plot[2] %in% tSNEObj$val@data@Dimnames[[1]]){
    for(i in features.plot)
    {  
      gene <- i
      # Create expression data frame for gene in long format
      df <- data.frame(Expression = tSNEObj$val@data[gene,], Gene = gene)
      # Add Cell column
      df$Cell = rownames(df)
      # Merge expression data frame to tSNE data frame
      df <- merge(df, tSNEmatrix$val, by = "Cell")

      write.table(df, "CoExpr.txt", sep = "\t", col.names = FALSE, row.names = FALSE, append = TRUE)
    }
    p = 1
    CoExprValue$val$Genes <- NULL
  } else if(features.plot[1] %in% tSNEObj$val@data@Dimnames[[1]]){
    noGenes[length(noGenes)+1] = features.plot[2]
    p = 0
  } else {
   noGenes[length(noGenes)+1] = features.plot[1]
   p = 0 
  }

GenesAbsent$val <- noGenes
if(p == 1){CoExprPLOT()}
