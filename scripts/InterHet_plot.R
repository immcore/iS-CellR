#!/usr/bin/env Rscript

if(!isS4(tSNEObj$val) || (mode$n == 0 && mode$m == 1) || "defLabels" %in% isolate(input$clustLabels) || input$changeLabels)
{
  if(!isS4(tSNEClusters$val) || (mode$n == 0 && mode$m == 1))
  {
    InterHetdata <- ProjectPCA(object = seuratObject$val, do.print = FALSE)
 
    InterHetdata <- FindClusters(object = InterHetdata, reduction.type = "pca", dims.use = 1:10, 
                            resolution = 0.6, print.output = 0, save.SNN = TRUE)

    # Make df.tsne global 
    tSNEClusters$val <- InterHetdata
  }

  if(!isS4(tSNEObj$val) || (mode$n == 0 && mode$m == 1))
  {
    InterHetdata <- RunTSNE(object = tSNEClusters$val, dims.use = 1:10, do.fast = TRUE)

    tSNEObj$val <- InterHetdata
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

if (file.exists("GeneSet1_avgExpr.txt")) file.remove("GeneSet1_avgExpr.txt")
if (file.exists("GeneSet2_avgExpr.txt")) file.remove("GeneSet2_avgExpr.txt")

Genes <- GeneListglob$val
GeneSet1.list <- as.vector(Genes[,1])
for(i in 1:length(GeneSet1.list)){
gene <- GeneSet1.list[i]
# Create expression data frame for gene in long format
df <- data.frame(Cell= names(tSNEObj$val@ident), Expression = tSNEObj$val@data[gene,], Gene = gene)
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

write.table(GeneSet1, "GeneSet1_avgExpr.txt", sep = "\t", col.names = FALSE, row.names = FALSE, append = TRUE)
}

GeneSet1_avgExpr <- read.delim("GeneSet1_avgExpr.txt", sep="\t", header=F)
names(GeneSet1_avgExpr) <- c("Cell", "Expression")
keys <- colnames(GeneSet1_avgExpr)[!grepl('Expression',colnames(GeneSet1_avgExpr))]
X <- as.data.table(GeneSet1_avgExpr)
GeneSet1 <- X[,list(GeneSet1avgExpr = mean(Expression)),keys]
##########################
GeneSet2.list <- as.vector(Genes[,2])
#length(GeneSet2.list)
for(i in 1:length(GeneSet2.list)){
  gene <- GeneSet2.list[i]
  # Create expression data frame for gene in long format
  df <- data.frame(Cell= names(tSNEObj$val@ident), Expression = tSNEObj$val@data[gene,], Gene = gene)
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
  
  write.table(GeneSet2, "GeneSet2_avgExpr.txt", sep = "\t", col.names = FALSE, row.names = FALSE, append = TRUE)
}

GeneSet2_avgExpr <- read.delim("GeneSet2_avgExpr.txt", sep="\t", header=F)
names(GeneSet2_avgExpr) <- c("Cell", "Expression")
keys <- colnames(GeneSet2_avgExpr)[!grepl('Expression',colnames(GeneSet2_avgExpr))]
X <- as.data.table(GeneSet2_avgExpr)
GeneSet2 <- X[,list(GeneSet2avgExpr= mean(Expression)),keys]

## Combine GeneSet1 and GeneSet2 average expression. Tirosh Fig3 A
GeneSet1_GeneSet2 <- merge(GeneSet1, GeneSet2, by = "Cell")

if (file.exists("GeneSet1_avgExpr.txt")) file.remove("GeneSet1_avgExpr.txt")
if (file.exists("GeneSet2_avgExpr.txt")) file.remove("GeneSet2_avgExpr.txt")

DownloadPlot$val$InterHet <- ggplot(GeneSet1_GeneSet2, aes(GeneSet1avgExpr, GeneSet2avgExpr)) + 
  geom_point(size=4, aes(GeneSet1avgExpr, GeneSet2avgExpr, colour = Cell)) + theme(legend.position="none") +
  xlab("Gene set1") + ylab("Gene set2") + 
  geom_text_repel(aes(GeneSet1avgExpr, GeneSet2avgExpr, label=Cell), size = 5) + ggtitle("Inter-sample heterogeneity")
