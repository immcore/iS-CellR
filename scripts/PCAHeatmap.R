#!/usr/bin/env Rscript

nCells <- input$PCAHeatCells
dims <- input$PCAHeatuse1

PCAheatdata <- ProjectPCA(object = seuratObject$val, do.print = FALSE)

PCHeatmap(object = PCAheatdata, pc.use = dims, cells.use = nCells, do.balanced = TRUE, label.columns = FALSE, do.return = TRUE)

DownloadPlot$val$PCAHeatmap <- recordPlot()