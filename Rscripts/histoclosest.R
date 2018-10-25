TMOhisto <- read.table("TMO_0.zerocounts", header = T, fill = T, row.names = NULL)
head(TMOhisto)
#par(mfrow=c(1,1))
hist(TMOhisto$NumberOfOccurance, xlab = "TMO", col="grey", breaks = 20, freq = NULL, main = "", plot = TRUE)
