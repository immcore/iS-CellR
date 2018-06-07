#!/usr/bin/env Rscript

expr = 0
countFiletmp <- as.data.frame(countFileQC$val)

freq <- as.numeric((rowSums(countFiletmp > expr)*100)/ncol(countFiletmp)) # Calculate frequency of cells in which each gene expressed
mean <- rowMeans(countFiletmp) # Calculate mean expression across all cells for each gene
df.freqmean <- cbind(freq, mean)

df.freqmean <- as.data.frame(df.freqmean)

genes <- dim(countFiletmp)[1]
cells <- dim(countFiletmp)[2]

DownloadPlot$val$Summaryplot  <- ggplot(data = df.freqmean, aes(x = mean,y = freq)) + geom_point(shape = 21, size = 1) + 
		geom_hline(yintercept = 50, colour = "red", linetype="dashed") + geom_smooth() +
		xlab(paste0("Mean expression (", genes ," genes)")) + 
		ylab(paste0("Frequency of expression (% of ", cells, " cells)")) + 
		theme(legend.title = element_blank()) +
           ggtitle(paste0("Frequency of expression vs Mean expression")) + theme(plot.title = element_text(color="black", size=14, face="bold"))

