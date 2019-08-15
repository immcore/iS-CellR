#!/usr/bin/env Rscript

nCells <- input$PCAHeatCells
dims <- input$PCAHeatuse1

#PCAheatdata <- ProjectPCA(object = scObject$val, do.print = FALSE)

#PCHeatmap(object = PCAheatdata, pc.use = dims, cells.use = nCells, do.balanced = TRUE, label.columns = FALSE, do.return = TRUE)

DownloadPlot$val$PCAHeatmap <- DimHeatmap(scObject$val, dims = dims, cells = nCells, balanced = TRUE)