#!/usr/bin/env Rscript

# Detection of variable genes across the single cells
scObject$val <- FindVariableFeatures(object = scObject$val, selection.method = 'mean.var.plot', mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(0.05, Inf))

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(scObject$val), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(scObject$val) + theme(legend.position="top")
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE) + theme(legend.position="top")
CombinePlots(plots = list(plot1, plot2))

