#!/usr/bin/env Rscript
if(isS4(scObject$val) || !"defLabels" %in% isolate(input$clustLabels) || input$changeLabels)
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

CoExprPLOT <- function(Dim1,Dim2,CoExpr_file){
    CoExpr <- CoExpr_file %>% dplyr::select(Cell, Expression, Gene, UMAP1, UMAP2, FItSNE_1, FItSNE_2, Cluster)

    CoExpr <- CoExpr_file
    names(CoExpr) <- c("Cell", "Expression", "Gene", "UMAP1", "UMAP2", "FItSNE_1", "FItSNE_2", "Cluster")
    #names(CoExpr) <- c("Cell", "Expression", "Gene", "UMAP1", "UMAP2", "Cluster")

if("UMAP" %in% isolate(dimPkg$val)) {
  CoExpr <- CoExpr[,c("Cell", "Expression", "Gene", "UMAP1", "UMAP2", "Cluster")]
} else {
  CoExpr <- CoExpr[,c("Cell", "Expression", "Gene", "FItSNE_1", "FItSNE_2", "Cluster")]
}

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

if(input$ExprCutoff > 0){
  Gene1 <- as.data.frame(subset(CoExpr, Gene == features.plot[1] & Expression >= input$ExprCutoff))
  Gene2 <- as.data.frame(subset(CoExpr, Gene == features.plot[2] & Expression >= input$ExprCutoff))
} else {
  Gene1 <- as.data.frame(subset(CoExpr, Gene == features.plot[1] & Expression > input$ExprCutoff))
  Gene2 <- as.data.frame(subset(CoExpr, Gene == features.plot[2] & Expression > input$ExprCutoff))
}

Gene12 <- merge(Gene1, Gene2, by = "Cell")
Gene12 <- Gene12[,1:9]
names(Gene12) <- c("Cell", "Expression.x", "Gene.x", Dim1, Dim2, "Cluster", "Celltype", "Expression.y", "Gene.y")

diff.Gene1 <- (!Gene1[,1] %in% Gene12[,1])
Gene1only <- Gene1[diff.Gene1,]
Gene1 <- Gene1only

diff.Gene2 <- (!Gene2[,1] %in% Gene12[,1])
Gene2only <- Gene2[diff.Gene2,]
Gene2 <- Gene2only

G1 <- paste0(features.plot[1]," Expr")
G2 <- paste0(features.plot[2]," Expr")
G12 <- paste0("Both Expr")

Gene12_cellCount <- function() {
  if("UMAP" %in% isolate(dimPkg$val)) {
  Gene12 %>%
    dplyr::group_by(Cluster) %>%
    dplyr::summarize(x = median(x = UMAP1), y = median(x = UMAP2)) -> labCluster
  } else {
  Gene12 %>%
    dplyr::group_by(Cluster) %>%
    dplyr::summarize(x = median(x = FItSNE_1), y = median(x = FItSNE_2)) -> labCluster
  }
  #ClustID <- CoExpr[tail(seq_along(CoExpr),2)] # Select last 2 columns
  ClustID <- Gene12 %>% dplyr::select(Cluster, Cell)
  #ClustID <- Gene12[tail(seq_along(Gene12),4)] # Select last 2 columns
  ClustID <- ClustID[,c(1,2)]
  lables <- as.data.frame(table(ClustID))
  lables <- lables[order(-lables$Freq), ] # count freq of clusters with celltype 
  ClustLab <- subset(lables, !duplicated(Cluster)) # remove duplicate or low counts
  labCount <- as.data.frame(table(ClustID[,1])) # count freq on cluster
  labCount <- labCount[order(-labCount$Freq), ]
  ClustLab <- cbind(ClustLab,labCount) 
  ClustInfo <- merge(labCluster, ClustLab, by = "Cluster")
  ClustInfo <- ClustInfo[,c(1,2,3,4,7)]
  names(ClustInfo) <- c("Cluster", Dim1, Dim2,  "Celltype", "Freq")
  Clustd <- ClustInfo[,c(1,5)]
  df.Gene12 <- merge(Gene12, Clustd, by = "Cluster")
  #write.table(df.Gene12, "features12.txt", sep = "\t", col.names = TRUE, row.names = FALSE)
  return(df.Gene12)
}

Gene1_cellCount <- function() {
  if("UMAP" %in% isolate(dimPkg$val)) {
  Gene1 %>%
    dplyr::group_by(Cluster) %>%
    dplyr::summarize(x = median(x = UMAP1), y = median(x = UMAP2)) -> labCluster
  } else {
  Gene1 %>%
    dplyr::group_by(Cluster) %>%
    dplyr::summarize(x = median(x = FItSNE_1), y = median(x = FItSNE_2)) -> labCluster
  }

  ClustID <- Gene1 %>% dplyr::select(Cluster, Cell)
  #ClustID <- Gene1[tail(seq_along(Gene1),2)] # Select last 2 columns
  lables <- as.data.frame(table(ClustID))
  lables <- lables[order(-lables$Freq), ] # count freq of clusters with celltype 
  
  ClustLab <- subset(lables, !duplicated(Cluster)) # remove duplicate or low counts
  labCount <- as.data.frame(table(ClustID[,1])) # count freq on cluster
  labCount <- labCount[order(-labCount$Freq), ]
  ClustLab <- cbind(ClustLab,labCount) 
  ClustInfo <- merge(labCluster, ClustLab, by = "Cluster")
  ClustInfo <- ClustInfo[,c(1,2,3,4,7)]
  names(ClustInfo) <- c("Cluster", Dim1, Dim2,  "Celltype", "Freq")
  Clustd <- ClustInfo[,c(1,5)]
  df.Gene1 <- merge(Gene1, Clustd, by = "Cluster")
  #write.table(df.Gene1, "features1.txt", sep = "\t", col.names = TRUE, row.names = FALSE)
  return(df.Gene1)
  #write.table(df.Gene1, "features1.txt", sep = "\t", col.names = TRUE, row.names = FALSE)
}

Gene2_cellCount <- function(){
  if("UMAP" %in% isolate(dimPkg$val)) {
  Gene2 %>%
    dplyr::group_by(Cluster) %>%
    dplyr::summarize(x = median(x = UMAP1), y = median(x = UMAP2)) -> labCluster
  } else {
  Gene2 %>%
    dplyr::group_by(Cluster) %>%
    dplyr::summarize(x = median(x = FItSNE_1), y = median(x = FItSNE_2)) -> labCluster
  }
  
  ClustID <- Gene2 %>% dplyr::select(Cluster, Cell)
  #ClustID <- Gene2[tail(seq_along(Gene2),2)] # Select last 2 columns
  lables <- as.data.frame(table(ClustID))
  lables <- lables[order(-lables$Freq), ] # count freq of clusters with celltype 
  ClustLab <- subset(lables, !duplicated(Cluster)) # remove duplicate or low counts
  labCount <- as.data.frame(table(ClustID[,1])) # count freq on cluster
  labCount <- labCount[order(-labCount$Freq), ]
  ClustLab <- cbind(ClustLab,labCount) 
  ClustInfo <- merge(labCluster, ClustLab, by = "Cluster")
  ClustInfo <- ClustInfo[,c(1,2,3,4,7)]
  names(ClustInfo) <- c("Cluster", Dim1, Dim2,  "Celltype", "Freq")
  Clustd <- ClustInfo[,c(1,5)]
  df.Gene2 <- merge(Gene2, Clustd, by = "Cluster")
  #write.table(df.Gene2, "features2.txt", sep = "\t", col.names = TRUE, row.names = FALSE)
  return(df.Gene2)
}


if("UMAP" %in% isolate(dimPkg$val)) {
  CoExpr %>%
    dplyr::group_by(Cluster) %>%
    dplyr::summarize(x = median(x = UMAP1), y = median(x = UMAP2)) -> labCluster
  } else {
  CoExpr %>%
    dplyr::group_by(Cluster) %>%
    dplyr::summarize(x = median(x = FItSNE_1), y = median(x = FItSNE_2)) -> labCluster
  }

  #ClustID <- CoExpr[tail(seq_along(CoExpr),2)] # Select last 2 columns
  ClustID <- CoExpr %>% dplyr::select(Cluster, Cell)
  #ClustID <- CoExpr[tail(seq_along(CoExpr),2)] # Select last 2 columns
  lables <- as.data.frame(table(ClustID))
  lables <- lables[order(-lables$Freq), ] # count freq of clusters with celltype 
  ClustLab <- subset(lables, !duplicated(Cluster)) # remove duplicate or low counts 
  labCount <- as.data.frame(table(ClustID[,1])) # count freq on cluster
  labCount <- labCount[order(-labCount$Freq), ]
  ClustLab <- cbind(ClustLab,labCount) 
  ClustInfo <- merge(labCluster, ClustLab, by = "Cluster")
  ClustInfo <- ClustInfo[,c(1,2,3,4,7)]
  names(ClustInfo) <- c("Cluster", Dim1, Dim2, "Celltype", "Freq")
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
  #ClustInfo$Cluster <- ClustInfo$Celltype
} 

  if(input$SwitchLabelCoExpr == "TRUE") {
        DownloadPlot$val$CoExprplot <- ggplot(data = CoExpr, aes_string(Dim1, Dim2)) + 
           geom_point(size=1, colour = "grey", aes(text=sprintf("Cell: %s<br>Cluster: %s [nCells: %s]<br>Gene: %s<br>Expr: %s", Cell, Cluster, Freq, Gene, Expression))) +
           geom_text(data=ClustInfo, mapping = aes(label = Cluster), size = 5, colour="black") +
           theme(legend.title = element_blank()) +
           ggtitle(paste0(features.plot[1]," x ",features.plot[2]," ","[Expr threshold: ",input$ExprCutoff,"]")) + 
           theme(plot.title = element_text(color="black", size=14, face="bold")) +
           labs(subtitle = paste0(features.plot[1],": ",signif(CoExprValue$val$Gene1Min, digits=3),"(min), ",signif(CoExprValue$val$Gene1Max,digits=3),"(max), ",signif(CoExprValue$val$Gene1Mean,digits=3),"(mean)\n",features.plot[2],": ",signif(CoExprValue$val$Gene2Min,digits=3),"(min), ",signif(CoExprValue$val$Gene2Max,digits=3),"(max), ",signif(CoExprValue$val$Gene2Mean,digits=3),"(mean)\n", "Total cells: ",nrow(Gene12)," (CoExpr), ", nrow(Gene1)," (",features.plot[1],"), ", nrow(Gene2)," (",features.plot[2],")\n")) +
  theme(plot.subtitle=element_text(size=11, hjust=0.5, color="black")) 

        DownloadPlot$val$CoExprplot1 <- ggplot(data = CoExpr, aes_string(Dim1, Dim2)) + 
           geom_point(size=1, colour = "grey", aes(text=sprintf("Cell: %s<br>Cluster: %s [nCells: %s]<br>Gene: %s<br>Expr: %s", Cell, Cluster, Freq, Gene, Expression))) +
           geom_text(data=ClustInfo, mapping = aes(label = Cluster), size = 5, colour="black") +
           theme(legend.title = element_blank()) +
           ggtitle(paste0(features.plot[1]," x ",features.plot[2]," ","[Expr threshold: ",input$ExprCutoff,"]")) + 
           theme(plot.title = element_text(color="black", size=14, face="bold")) + 
           labs(subtitle = paste0(features.plot[1],": ",signif(CoExprValue$val$Gene1Min, digits=3),"(min), ",signif(CoExprValue$val$Gene1Max,digits=3),"(max), ",signif(CoExprValue$val$Gene1Mean,digits=3),"(mean)\n",features.plot[2],": ",signif(CoExprValue$val$Gene2Min,digits=3),"(min), ",signif(CoExprValue$val$Gene2Max,digits=3),"(max), ",signif(CoExprValue$val$Gene2Mean,digits=3),"(mean)\n", "Total cells: ",nrow(Gene12)," (CoExpr), ", nrow(Gene1)," (",features.plot[1],"), ", nrow(Gene2)," (",features.plot[2],")\n")) +
  theme(plot.subtitle=element_text(size=11, hjust=0.5, color="black")) 
    } else {
        DownloadPlot$val$CoExprplot <- ggplot(data = CoExpr, aes_string(Dim1, Dim2)) + 
           geom_point(size=1, colour = "grey", aes(text=sprintf("Cell: %s<br>Cluster: %s [nCells: %s]<br>Gene: %s<br>Expr: %s", Cell, Cluster, Freq, Gene, Expression))) +
           theme(legend.title = element_blank()) +
           ggtitle(paste0(features.plot[1]," x ",features.plot[2]," ","[Expr threshold: ",input$ExprCutoff,"]")) + 
           theme(plot.title = element_text(color="black", size=14, face="bold")) +
           labs(subtitle = paste0(features.plot[1],": ",signif(CoExprValue$val$Gene1Min, digits=3),"(min), ",signif(CoExprValue$val$Gene1Max,digits=3),"(max), ",signif(CoExprValue$val$Gene1Mean,digits=3),"(mean)\n",features.plot[2],": ",signif(CoExprValue$val$Gene2Min,digits=3),"(min), ",signif(CoExprValue$val$Gene2Max,digits=3),"(max), ",signif(CoExprValue$val$Gene2Mean,digits=3),"(mean)\n", "Total cells: ",nrow(Gene12)," (CoExpr), ", nrow(Gene1)," (",features.plot[1],"), ", nrow(Gene2)," (",features.plot[2],")\n")) +
  theme(plot.subtitle=element_text(size=11, hjust=0.5, color="black"))

        DownloadPlot$val$CoExprplot1 <- ggplot(data = CoExpr, aes_string(Dim1, Dim2)) + 
           geom_point(size=1, colour = "grey", aes(text=sprintf("Cell: %s<br>Cluster: %s [nCells: %s]<br>Gene: %s<br>Expr: %s", Cell, Cluster, Freq, Gene, Expression))) +
           theme(legend.title = element_blank()) +
           ggtitle(paste0(features.plot[1]," x ",features.plot[2]," ","[Expr threshold: ",input$ExprCutoff,"]")) + 
           theme(plot.title = element_text(color="black", size=14, face="bold")) + 
           labs(subtitle = paste0(features.plot[1],": ",signif(CoExprValue$val$Gene1Min, digits=3),"(min), ",signif(CoExprValue$val$Gene1Max,digits=3),"(max), ",signif(CoExprValue$val$Gene1Mean,digits=3),"(mean)\n",features.plot[2],": ",signif(CoExprValue$val$Gene2Min,digits=3),"(min), ",signif(CoExprValue$val$Gene2Max,digits=3),"(max), ",signif(CoExprValue$val$Gene2Mean,digits=3),"(mean)\n", "Total cells: ",nrow(Gene12)," (CoExpr), ", nrow(Gene1)," (",features.plot[1],"), ", nrow(Gene2)," (",features.plot[2],")\n")) +
  theme(plot.subtitle=element_text(size=11, hjust=0.5, color="black"))
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
  
  #if (file.exists(CoExpr_file)) file.remove(CoExpr_file)
}

Dim1 <- paste(dimPkg$val,"1",sep="")
Dim2 <- paste(dimPkg$val,"2",sep="")

# CoExpr gene of interest
features.plot <- as.character(unlist(strsplit(input$CoExprGenes," "))) 

DownloadPlot$val$CoExprplot <- NULL
CoExprValue$val$Gene1Ab <- NULL
CoExprValue$val$Gene2Ab <- NULL
CoExpr_file <- data.frame()
#now <- Sys.time()
#CoExpr_file <- paste0("/tmpfiles/",str_replace_all(input$username,"\\s+","_"),"_",format(now,"%Y%m%d_"),"CoExpr.txt")

#if (file.exists(CoExpr_file)) file.remove(CoExpr_file)

noGenes <- c()
p = 0

if(features.plot[1] %in% GetAssayData(object = scObject$val)@Dimnames[[1]] && features.plot[2] %in% GetAssayData(object = scObject$val)@Dimnames[[1]]){
    for(i in features.plot)
    {  
      gene <- i
      # Create expression data frame for gene in long format
      df <- data.frame(Expression = GetAssayData(object = scObject$val)[gene,], Gene = gene)
      # Add Cell column
      df$Cell = rownames(df)
      # Merge expression data frame to tSNE data frame
      df <- merge(df, tSNEmatrix$val, by = "Cell")

      CoExpr_file <- rbind(CoExpr_file, df)

      #write.table(df, CoExpr_file, sep = "\t", col.names = FALSE, row.names = FALSE, append = TRUE)
    }
    p = 1
    CoExprValue$val$Genes <- NULL
  } else if(features.plot[1] %in% GetAssayData(object = scObject$val)@Dimnames[[1]]){
    noGenes[length(noGenes)+1] = features.plot[2]
    p = 0
  } else {
   noGenes[length(noGenes)+1] = features.plot[1]
   p = 0 
  }

GenesAbsent$val <- noGenes
if(p == 1){CoExprPLOT(Dim1,Dim2,CoExpr_file)}
