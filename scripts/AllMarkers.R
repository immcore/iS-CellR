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

if("useheader" %in% isolate(input$clustLabels)) {
  current.clustID <- as.data.frame(scObject$val@ident)
  setDT(current.clustID, keep.rownames = TRUE)[]
  colnames(current.clustID) <- c("clust", "ident")
  current.clustID$clust <- gsub("\\..*|_.*|-.*", "", current.clustID$clust)
  current.clustID %>% group_by(ident,clust) %>% tally() -> lables
  new.lables <- as.data.frame(lables)
  new.lables <- new.lables[order(-new.lables$n), ]
  ClustLab <- subset(new.lables, !duplicated(ident))
  new.ident <- as.vector(ClustLab$ident)
  new.clust <- as.vector(ClustLab$clust)

  scObject$val@ident <- plyr::mapvalues(x = scObject$val@ident, from = new.ident, to = new.clust)
}

if(is.null(scObjAllmarkers$val)){
    # find markers for every cluster compared to all remaining cells, report only the positive ones
    scObjAllmarkers$val <- FindAllMarkers(object = scObject$val, test.use = "MAST")
}

if(input$Topmarkers > 0 ){
    Allmarkers <- scObjAllmarkers$val %>% group_by(cluster) %>% top_n(input$Topmarkers, avg_logFC)
    Allmarkers$val = Allmarkers
} else {
    scObjAllmarkers$val
} 
