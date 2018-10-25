#!/usr/bin/env Rscript

# Detection of variable genes across the single cells
VarGenes <- FindVariableGenes(object = scObject$val, mean.function = ExpMean, dispersion.function = LogVMR, 
                          x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

