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

#now <- Sys.time()
#GeneSet1_file <- paste0("/tmpfiles/",str_replace_all(input$username,"\\s+","_"),"_",format(now,"%Y%m%d_"),"GeneSet1_avgExpr.txt")
#GeneSet2_file <- paste0("/tmpfiles/",str_replace_all(input$username,"\\s+","_"),"_",format(now,"%Y%m%d_"),"GeneSet2_avgExpr.txt")
GeneSet1_file <- data.frame()
GeneSet2_file <- data.frame()
#if (file.exists(GeneSet1_file)) file.remove(GeneSet1_file)
#if (file.exists(GeneSet2_file)) file.remove(GeneSet2_file)

Genes <- GeneListglob$val
GeneSet1.list <- as.vector(Genes[,1])
for(i in 1:length(GeneSet1.list)){
  gene <- GeneSet1.list[i]

  if(gene %in% GetAssayData(object = scObject$val)@Dimnames[[1]]){
    # Create expression data frame for gene in long format
    df <- data.frame(Cell= names(Idents(object = scObject$val)), Expression = GetAssayData(object = scObject$val)[i,], Gene = gene)
    # Merge expression data frame to tSNE data frame
    df <- merge(df, tSNEmatrix$val, by = "Cell")
    # Extract sample name from Cell and add column

    df$Cell <- gsub("\\..*|_.*|-.*", "", df$Cell)
    if("useheader" %in% isolate(input$clustLabels)) {
      GeneSet1_gene <- df[c(1,2)]
    } else {
      GeneSet1_gene <- df[c(6,2)]
    }
    ## count average expression of gene for sample
    keys <- colnames(GeneSet1_gene)[!grepl('Expression',colnames(GeneSet1_gene))]
    X <- as.data.table(GeneSet1_gene)
    GeneSet1 <- X[,list(mm= mean(Expression)),keys]

    GeneSet1_file <- rbind(GeneSet1_file, GeneSet1)
    #write.table(GeneSet1, GeneSet1_file, sep = "\t", col.names = FALSE, row.names = FALSE, append = TRUE)
  }
}

#GeneSet1_avgExpr <- read.delim(GeneSet1_file, sep="\t", header=F)
GeneSet1_avgExpr <- GeneSet1_file
names(GeneSet1_avgExpr) <- c("Cell", "Expression")
keys <- colnames(GeneSet1_avgExpr)[!grepl('Expression',colnames(GeneSet1_avgExpr))]
X <- as.data.table(GeneSet1_avgExpr)
GeneSet1 <- X[,list(GeneSet1avgExpr = mean(Expression)),keys]
##########################
GeneSet2.list <- as.vector(Genes[,2])
#length(GeneSet2.list)
for(i in 1:length(GeneSet2.list)){
  gene <- GeneSet2.list[i]
  
  if(gene %in% GetAssayData(object = scObject$val)@Dimnames[[1]]){
    # Create expression data frame for gene in long format
    df <- data.frame(Cell= names(Idents(object = scObject$val)), Expression = GetAssayData(object = scObject$val)[i,], Gene = gene)
    # Merge expression data frame to tSNE data frame
    df <- merge(df, tSNEmatrix$val, by = "Cell")
    # Extract sample name from Cell and add column
    #df$sample <- sapply(strsplit(df$Cell,"_"), `[`, 1)
    df$Cell <- gsub("\\..*|_.*|-.*", "", df$Cell)
    
    if("useheader" %in% isolate(input$clustLabels)) {
      GeneSet2_gene <- df[c(1,2)]
      } else {
      GeneSet2_gene <- df[c(6,2)]  
      }
     ## count average expression of gene for sample
    keys <- colnames(GeneSet2_gene)[!grepl('Expression',colnames(GeneSet2_gene))]
    X <- as.data.table(GeneSet2_gene)
    GeneSet2 <- X[,list(mm= mean(Expression)),keys]
    
    GeneSet2_file <- rbind(GeneSet2_file, GeneSet2)
    #write.table(GeneSet2, GeneSet2_file, sep = "\t", col.names = FALSE, row.names = FALSE, append = TRUE)
  }
}

#GeneSet2_avgExpr <- read.delim(GeneSet2_file, sep="\t", header=F)
GeneSet2_avgExpr <- GeneSet2_file
names(GeneSet2_avgExpr) <- c("Cell", "Expression")
keys <- colnames(GeneSet2_avgExpr)[!grepl('Expression',colnames(GeneSet2_avgExpr))]
X <- as.data.table(GeneSet2_avgExpr)
GeneSet2 <- X[,list(GeneSet2avgExpr= mean(Expression)),keys]

## Combine GeneSet1 and GeneSet2 average expression. Tirosh Fig3 A
GeneSet1_GeneSet2 <- merge(GeneSet1, GeneSet2, by = "Cell")

#if (file.exists(GeneSet1_file)) file.remove(GeneSet1_file)
#if (file.exists(GeneSet2_file)) file.remove(GeneSet2_file)
GeneSet1_file <- data.frame()
GeneSet2_file <- data.frame()

DownloadPlot$val$InterHet <- ggplot(GeneSet1_GeneSet2, aes(GeneSet1avgExpr, GeneSet2avgExpr)) + 
  geom_point(size=4, aes(GeneSet1avgExpr, GeneSet2avgExpr, colour = Cell)) + theme(legend.position="none") +
  xlab("Gene set1") + ylab("Gene set2") + 
  
  geom_text_repel(aes(GeneSet1avgExpr, GeneSet2avgExpr, label=Cell), size = 5) + ggtitle("Inter-sample heterogeneity")
