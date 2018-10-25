#!/usr/bin/env Rscript

#if(!is.null(DistinguishMarkers$val))
#{
  df.volcano <- as.data.frame(DistinguishMarkers$val)
  fwrite(df.volcano, "up.txt", sep = "\t", col.names = TRUE, row.names = TRUE)

  df.volcano$Gene <- rownames(df.volcano)
  dim(df.volcano)
  # add a grouping column; default value is "not significant"
 ############# For GO annotation - Significantly UP in ident.1 pval < 0.05 and logFC > 1
#df.Expr$Gene <- rownames(df.Expr)
up.ident1 <- subset(df.volcano, df.volcano$p_val_adj<0.05 & df.volcano$avg_logFC>1)

############# For GO annotation - Significantly Down in ident.1 pval < 0.05 and logFC < 1
#df.Expr$Gene <- rownames(df.Expr)
#down.ident1 <- subset(df.volcano, p_val_adj<0.05 & avg_logFC<(-1))


DownloadPlot$val$Upplot <- ggplotly(ggplot(data = up.ident1, aes(x= pct.1, y = pct.2)) + 
  geom_point(color = "blue", size = 3, mapping = aes(text=sprintf("Gene: %s", Gene))) + ggtitle("Up regulated genes") +
  xlab("CT83pos cells") + ylab("CT83neg cells") + 
  scale_x_continuous(limits = c(0,1.00), breaks = c(0,0.25,0.50,0.75,1.00),
                     labels = c("0%","25%","50%","75%","100%")) +
  scale_y_continuous(limits = c(0,1.00), breaks = c(0,0.25,0.50,0.75,1.00),
                     labels = c("0%","25%","50%","75%","100%")),
  tooltip = c("text"))
