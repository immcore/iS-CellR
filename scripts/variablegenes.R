#!/usr/bin/env Rscript

# Detection of variable genes across the single cells

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(scObject$val), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(scObject$val)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))

