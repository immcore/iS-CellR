#!/usr/bin/env Rscript

if(isS4(scObject$val) || "defLabels" %in% isolate(input$clustLabels) || input$changeLabels)
{
	################## for Custom labels ################
  if("customLabels" %in% isolate(input$clustLabels)) {
    cluster.ids <- as.character(unlist(ClusterLabInfo$val[,1]))#, decreasing = FALSE)#c("CD4", "Bcells", "CD8cells",
    new.cluster.ids <- as.character(unlist(ClusterLabInfo$val[,2]))#c("CD4", "Bcells", "CD8cells", 
    if(input$changeLabels){
      scObject$val@ident <- plyr::mapvalues(x = scObject$val@ident, from = cluster.ids, to = new.cluster.ids)
    }
    CInfo <- cbind(cluster.ids,new.cluster.ids)
    ClusterLabInfo$val <- CInfo
  }

  if("defLabels" %in% isolate(input$clustLabels)) {
    if(!is.null(dfcluster.ids$val)){
      new.cluster.ids <- as.character(unlist(ClusterLabInfo$val[,2]))
      current.ids <- as.character(unlist(ClusterLabInfo$val[,1]))#, decreasing = FALSE)#c("CD4", "Bcells", "CD8cells",
      scObject$val@ident <- plyr::mapvalues(x = scObject$val@ident, from = new.cluster.ids, to = current.ids)
    } else {
      new.cluster.ids = ""
      current.ids <- sort(as.character(unique(scObject$val@ident)), decreasing = FALSE)
    }
    cluster.ids <- current.ids
    CInfo <- cbind(cluster.ids,new.cluster.ids)
    ClusterLabInfo$val <- CInfo 
    dfcluster.ids$val <- cluster.ids      
  }

mode$m <- 0
}

if(input$SwitchHeatmap == "TRUE"){
    if(is.null(scObjAllmarkers$val)){
      # find markers for every cluster compared to all remaining cells, report only the positive ones
        scObjAllmarkers$val <- FindAllMarkers(object = scObject$val, only.pos = TRUE, min.pct = 0.25, 
                                              thresh.use = 0.25)
    }
    topN <- scObjAllmarkers$val %>% group_by(cluster) %>% top_n(input$TopDiffHeatmap, avg_logFC)
    
} else {


      # find all markers of cluster 1
      Markers_Clust1 <- FindMarkers(object = scObject$val, ident.1 = input$ClusterSel, min.pct = 0.25)
      #print(x = head(x = cluster1.markers, n = 5))

      # find all markers distinguishing cluster 5 from clusters 0 and 3
      Distin_Markers <- FindMarkers(object = scObject$val, ident.1 = input$ClusterSel, ident.2 = c(input$ClusterSel, input$ClusterSel), min.pct = 0.25)
      #print(x = head(x = cluster5.markers, n = 5))



    if(!is.null(input$HeatmapGenes)){
        #Five visualizations of marker gene expression
        features.plot <- as.character(unlist(strsplit(input$HeatmapGenes,","))) 
        #Single cell heatmap of gene expression

        noGenes <- c()
        Genes <- c()
    
        for (i in features.plot){
          if(i %in% rownames(scObject$val@scale.data))
          {
              Genes[length(Genes)+1] = i
          } else {
              noGenes[length(noGenes)+1] = i
          }
        }

        GenesAbsent$val <- noGenes

        if(length(GenesAbsent$val) == length(features.plot)) {
              plot.new()
              DownloadPlot$val$Heatmap <- NULL
        } else {
              DownloadPlot$val$Heatmap <- DoHeatmap(object = scObject$val, genes.use = Genes, 
                                            slim.col.label = TRUE, remove.key = TRUE)
        }
    } else {
        plot.new()
        DownloadPlot$val$Heatmap <- NULL
    }
}