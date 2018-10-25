#!/usr/bin/env Rscript

if(is.null(DownloadPlot$val$Downplot)){

  df.volcano <- as.data.frame(DistinguishMarkers$val)

  df.volcano$Gene <- rownames(df.volcano)

  # add a grouping column; default value is "not significant"
  ############# For GO annotation - Significantly Down in ident.1 pval < 0.05 and logFC < 1
  #df.Expr$Gene <- rownames(df.Expr)
  down.ident1 <- subset(df.volcano, p_val_adj<0.05 & avg_logFC<(-1))

  if(!is.null(down.ident1)){
    DownloadPlot$val$Downplot <- ggplotly(ggplot(data = down.ident1, aes(x= pct.1, y = pct.2, colour = pct.1)) + 
      geom_point(size = 2, mapping = aes(text=sprintf("Gene: %s", Gene))) + ggtitle(paste0("Down regulated genes in Cluster: ",input$MarkersCluster1)) +
      xlab("Group 1 cells") + ylab("Group 2 cells") + 
      scale_colour_viridis(option = "plasma", direction = -1) + 
      theme(legend.position = "none") +
      scale_x_continuous(limits = c(0,1.00), breaks = c(0,0.25,0.50,0.75,1.00),
                         labels = c("0%","25%","50%","75%","100%")) +
      scale_y_continuous(limits = c(0,1.00), breaks = c(0,0.25,0.50,0.75,1.00),
                         labels = c("0%","25%","50%","75%","100%")),
      tooltip = c("text"))
    } else {
      DownloadPlot$val$Downplot <- NULL
    }
} else {
  DownloadPlot$val$Downplot <- NULL
}