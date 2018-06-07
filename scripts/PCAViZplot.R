#!/usr/bin/env Rscript

nGenes <- input$PCAvizGenes
font <- input$PCAvizfont
plot_list <- list()
dims = input$VizPCAuse1:input$VizPCAuse2

df.pcs <- as.data.frame(PCAClustGlob$val@dr$pca@gene.loadings)

for (i in dims){
 
 	dataset <- df.pcs[DimTopGenes(
      object = PCAClustGlob$val,
      dim.use = i,
      reduction.type = "pca",
      num.genes = nGenes,
      use.full = FALSE,
      do.balanced = FALSE
    ),]

    plot_single <- ggplot(data = dataset, mapping = aes(dataset[,i], y=factor(1:nrow(x=dataset)))) +
			geom_point(size = 1, color = "red") +
			xlab(paste0("PC", i)) + ylab("Genes") +
			scale_y_discrete(labels = rownames(dataset)) +
			theme(axis.text = element_text(size=font)) +
			theme(axis.title.x = element_text(hjust=0.5)) +
			theme(axis.title.y = element_text(hjust=0.5)) +
			ggtitle(paste0("PC", i)) + theme(plot.title = element_text(color="black", size=14, face="bold")) + theme(legend.position="none") + theme(legend.title = element_blank())

	plot_list[[length(plot_list)+1]] <- plot_single

}

grid.arrange(grobs=plot_list, ncol=3)

DownloadPlot$val$Vizplot <- plot_grid(plotlist = plot_list, ncol = 3)

