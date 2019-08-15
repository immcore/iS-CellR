#!/usr/bin/env Rscript

#now <- Sys.time()
#GeneSet1_sCell_file <- paste0("/tmpfiles/",str_replace_all(input$username,"\\s+","_"),"_",format(now,"%Y%m%d_"),"GeneSet1_avgExpr_sCell.txt")
#GeneSet2_sCell_file <- paste0("/tmpfiles/",str_replace_all(input$username,"\\s+","_"),"_",format(now,"%Y%m%d_"),"GeneSet2_avgExpr_sCell.txt")
GeneSet1_sCell_file <- data.frame()
GeneSet2_sCell_file <- data.frame()
#if (file.exists(GeneSet1_sCell_file)) file.remove(GeneSet1_sCell_file)
#if (file.exists(GeneSet2_sCell_file)) file.remove(GeneSet2_sCell_file)

Genes <- GeneListglob$val

GeneSet1.list <- as.vector(Genes[,1])
length(GeneSet1.list)
for(i in 1:length(GeneSet1.list)){
  gene <- GeneSet1.list[i]
  # Create expression data frame for gene in long format
  if(gene %in% GetAssayData(object = scObject$val)@Dimnames[[1]]){
    df <- data.frame(Cell= names(Idents(object = scObject$val)), Expression = GetAssayData(object = scObject$val)[i,], Gene = gene)
    # Merge expression data frame to tSNE data frame
    df <- merge(df, tSNEmatrix$val, by = "Cell")
    
    GeneSet1_gene <- df[c(1,6,2)]

    GeneSet1_sCell_file <- rbind(GeneSet1_sCell_file, GeneSet1_gene)
    #write.table(GeneSet1_gene, GeneSet1_sCell_file, sep = "\t", col.names = FALSE, row.names = FALSE, append = TRUE)
  }
}

#GeneSet1_avgExprSC <- read.delim(GeneSet1_sCell_file, sep="\t", header=F)
GeneSet1_avgExprSC <- GeneSet1_sCell_file
names(GeneSet1_avgExprSC) <- c("Cell", "Cluster", "Expression")
keys <- colnames(GeneSet1_avgExprSC)[!grepl('Expression',colnames(GeneSet1_avgExprSC))]
X <- as.data.table(GeneSet1_avgExprSC)
GeneSet1SC <- X[,list(GeneSet1avgExprSC= mean(Expression)),keys]
##########################

GeneSet2.list <- as.vector(Genes[,2])
length(GeneSet2.list)
for(i in 1:length(GeneSet2.list)){
  gene <- GeneSet2.list[i]
  
  if(gene %in% GetAssayData(object = scObject$val)@Dimnames[[1]]){
    # Create expression data frame for gene in long format
    df <- data.frame(Cell= names(Idents(object = scObject$val)), Expression = GetAssayData(object = scObject$val)[i,], Gene = gene)
    # Merge expression data frame to tSNE data frame
    df <- merge(df, tSNEmatrix$val, by = "Cell")
    # Extract sample name from Cell and add column
    
    GeneSet2_gene <- df[c(1,6,2)]  
    
    GeneSet2_sCell_file <- rbind(GeneSet2_sCell_file, GeneSet2_gene)
    #write.table(GeneSet2_gene, GeneSet2_sCell_file, sep = "\t", col.names = FALSE, row.names = FALSE, append = TRUE)
  }
}

#GeneSet2_avgExprSC <- read.delim(GeneSet2_sCell_file, sep="\t", header=F)
GeneSet2_avgExprSC <- GeneSet2_sCell_file
names(GeneSet2_avgExprSC) <- c("Cell", "Cluster", "Expression")
keys <- colnames(GeneSet2_avgExprSC)[!grepl('Expression',colnames(GeneSet2_avgExprSC))]
X <- as.data.table(GeneSet2_avgExprSC)
GeneSet2SC <- X[,list(GeneSet2avgExprSC= mean(Expression)),keys]

GeneSet1_GeneSet2_SC <- merge(GeneSet1SC, GeneSet2SC, by = "Cell")

if("useheader" %in% isolate(input$clustLabels)) {
  GeneSet1_GeneSet2_SC$CellName <- GeneSet1_GeneSet2_SC$Cell
  GeneSet1_GeneSet2_SC$CellName <- gsub("\\..*|_.*|-.*", "", GeneSet1_GeneSet2_SC$CellName)
} else {
  GeneSet1_GeneSet2_SC$CellName <- GeneSet1_GeneSet2_SC[,2]
}

#if (file.exists(GeneSet1_sCell_file)) file.remove(GeneSet1_sCell_file)
#if (file.exists(GeneSet2_sCell_file)) file.remove(GeneSet2_sCell_file)
GeneSet1_sCell_file <- data.frame()
GeneSet2_sCell_file <- data.frame()

list.Celltypes <- as.character(unlist(strsplit(input$IntraHetGenes,","))) 

f1 <- list(
  family = "Arial, sans-serif",
  size = 14,
  color = "lightgrey"
)

xaxis <- list(
    title = "Gene list1",
    titlefont = f1,
    range = c(0, 1.25),
    showgrid = FALSE,
    zeroline = FALSE,
    showticklabels = TRUE,
    showline = TRUE
  )
yaxis <- list(
    title = "Gene list2",
    titlefont = f1,
    range = c(0, 1),
    showgrid = FALSE,
    zeroline = FALSE,
    showticklabels = TRUE,
    showline = TRUE
  )

m <- list(
  l = 50,
  r = 50,
  b = 50,
  t = 50,
  pad = 0
)

plot_list <- list()

plot_df <- data.frame()
noGenes <- c()

for (i in list.Celltypes) {
  if(i %in% GeneSet1_GeneSet2_SC$CellName)
  {
    plot_df <- rbind(plot_df, GeneSet1_GeneSet2_SC[GeneSet1_GeneSet2_SC$CellName == i,])
   } else{
    noGenes[length(noGenes)+1] = i
  }
}

GenesAbsent$val <- noGenes

Set1Max <- max(plot_df$GeneSet1avgExprSC)
Set2Max <- max(plot_df$GeneSet2avgExprSC)

if(Set1Max > Set2Max){
  xMax <- Set1Max
  yMax <- Set1Max
} else{
  xMax <- Set2Max
  yMax <- Set2Max
}

DownloadPlot$val$IntraHetplotly <- ggplotly(ggplot(data = plot_df, aes(GeneSet1avgExprSC, GeneSet2avgExprSC)) + 
      geom_point(mapping = aes(text=sprintf("Cell: %s<br>Cluster: %s", Cell, CellName), color = GeneSet2avgExprSC)) +
      ggtitle("Intra-sample heterogeneity") + theme(plot.title = element_text(color="black", size=14, face="bold")) + 
      scale_x_continuous(limits = c(0, xMax)) +
      scale_y_continuous(limits = c(0, yMax)) + facet_wrap(~CellName, scales = "free") + 
      theme(strip.background = element_blank(), strip.placement = "outside") +
      theme(strip.text = element_text(size = rel(1))) +
      xlab("Gene set1") + ylab("Gene set2") +
      theme(axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 2, 0, 0), "mm"), angle = 90)) +
      scale_color_gradientn(name="AvgExpr", colors= c("blue", "orange")) + theme(legend.position="top") + 
      labs(color='AvgExpr') +
      theme(plot.margin = unit(c(1,1,1,1), "cm"), panel.margin = unit(0, "cm")) +
      #theme(panel.spacing=unit(4, "cm")) +
      theme(legend.title = element_blank()),
      tooltip = c("text"))


DownloadPlot$val$IntraHetplot <- ggplot(data = plot_df, aes(GeneSet1avgExprSC, GeneSet2avgExprSC)) + 
      geom_point(mapping = aes(text=sprintf("Cell: %s<br>Cluster: %s", Cell, CellName), color = GeneSet2avgExprSC)) +
      ggtitle("Intra-sample heterogeneity") + theme(plot.title = element_text(color="black", size=14, face="bold")) + theme(legend.position="right") + 
      scale_x_continuous(limits = c(0, xMax)) +
      scale_y_continuous(limits = c(0, yMax)) + facet_wrap(~CellName, scales = "free") + 
      theme(strip.background = element_blank(), strip.placement = "outside") +
      theme(strip.text = element_text(size = rel(1))) +
      xlab("Gene set1") + ylab("Gene set2") +
      theme(axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 2, 0, 0), "mm"), angle = 90)) +
      scale_color_gradientn(colors= c("blue", "orange")) +
      theme(plot.margin = unit(c(1,1,1,1), "cm"), panel.margin = unit(0, "cm")) +
      #theme(panel.spacing=unit(4, "cm")) +
      theme(legend.title = element_blank())

