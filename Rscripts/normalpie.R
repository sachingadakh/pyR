library(plotly)

pie(table(TMO_0$Frequency), radius = 1.5)
pie(read.table("TMO_0.zerocounts"), radius = 1.5)
pie(TMO_0$Frequency, as.graphicsAnnot(TMO_0$Frequency),col = rainbow(6), main = "Frequency of repeats in TMO",radius = 1)
legend("topright", legend = TMO_0$NameOfRepeats, fill = rainbow(6), cex = 0.8)

TMOpiedata <- data.frame(TMO_0)
TMOpiedata 
lbls <- paste(names(TMOpiedata), TMOpiedata, sep = "")
lbls
pie(TMOpiedata, labels = lbls, col = rainbow(6) ,main = "Frequency of repeats in TMO", radius = 1)


df <- TMO_0
df$Frequency <- sapply(rownames(TMO_0), function(x) strsplit(x, split = "")[[1]][1])
plot.df <- df 
group_by(Frequency)   
summarize(count = n())

