#!/usr/bin/env Rscript

nGenes <- input$PCAvizGenes
font <- input$PCAvizfont
plot_list <- list()
dims = input$VizPCAuse1:input$VizPCAuse2

DownloadPlot$val$Vizplot <- VizDimLoadings(scObject$val, dims = dims, reduction = "pca", nfeatures = nGenes, col = "blue")


