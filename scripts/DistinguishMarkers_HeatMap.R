#!/usr/bin/env Rscript

if(is.null(DownloadPlot$val$DistHeatmap)){

  df.volcano <- as.data.frame(DistinguishMarkers$val)

  df.volcano$Gene <- rownames(df.volcano)

  genes <- subset(df.volcano, p_val_adj<0.05 & abs(avg_logFC)>1) # gives all up and down genes for corresponding cutoff 
  #genes <- subset(de_genes, p_val_adj<0.05 & avg_logFC>1) # gives all up genes in ident1 
  #genes <- subset(de_genes, p_val_adj<0.05 & avg_logFC<(-1)) # gives all down genes in ident1

  MarkersCluster2 <- as.character(unlist(strsplit(input$MarkersCluster2," "))) 
  MarkersCluster1 <- as.character(input$MarkersCluster1)
  cells <- c(MarkersCluster1, MarkersCluster2)

  Cells.use <- SubsetData(object = scObject$val, ident.use = cells)

  DownloadPlot$val$DistHeatmap <- DoHeatmap(object = Cells.use, genes.use = rownames(genes), rotate.key=TRUE, cex.row = 12,
                           slim.col.label = TRUE, remove.key = FALSE)

  DownloadPlot$val$DistHeatmap <- ggplotly(DownloadPlot$val$DistHeatmap + scale_fill_viridis(option = "viridis")) 

} else {
  DownloadPlot$val$DistHeatmap
}