genecount <- read.table("TMOannotated.TSV", header = T, fill = T, row.names = NULL)
par(mfrow= c(1,1))
gc1 <- table(genecount$geneId)
gc1
gc2 <- as.data.frame(table(genecount$geneId)) 
gc2
hist(gc2$Freq, xlab = "TMO", col = "grey", main="", breaks = 10)
hist(gc2$Freq, main = "TMO", col = "blue", breaks = c(0:10), xlab="Number of DB sites (bins)", ylab ="Gene frequency", ylim = c(0,100))
write.table(as.data.frame(gc2), file = "TMOcount", quote = FALSE, sep = "\t", row.names = F)

genecounttmy <- read.table("TMYannotated.TSV", header = T, fill = T, row.names = NULL)
par(mfrow= c(1,1))
gctmy <- table(genecounttmy$geneId)
gctmy
gctmy2 <- as.data.frame(table(genecounttmy$geneId))
gctmy2
hist(gctmy2$Freq, xlab = "TMY", col = "grey", main="", breaks = 10)
hist(gctmy2$Freq, main = "TMY", col = "red", breaks = c(0:20), xlab="Number of DB sites (bins)", ylab ="Gene frequency", ylim = c(0,600))
write.table(as.data.frame(gctmy2), file = "TMYcount", quote = FALSE, sep = "\t", row.names = F)

genecounttfo <- read.table("TFOannotated.TSV", header = T, fill = T, row.names = NULL)
par(mfrow= c(1,1))
gctfo <- table(genecounttfo$geneId)
gctfo
gctfo2 <- as.data.frame(table(genecounttfo$geneId))
gctfo2
hist(gctfo2$Freq, xlab = "TFO", col = "grey", main = "", breaks = 10)
hist(gctfo2$Freq, main = "TFO", col = "blue", breaks = c(0:10), xlab="Number of DB sites (bins)", ylab ="Gene frequency", ylim = c(0,200))
write.table(as.data.frame(gctfo2), file = "TFOcount", quote = FALSE, sep = "\t", row.names = F)

genecounttfy <- read.table("TFYannotated.TSV", header = T, fill = T, row.names = NULL)
par(mfrow= c(1,1))
gctfy <- table(genecounttfy$geneId)
gctfy
gctfy2 <- as.data.frame(table(genecounttfy$geneId))
gctfy2
hist(gctfy2$Freq, xlab = "TFY", col = "grey", main = "", breaks = 10)
hist(gctfy2$Freq, main = "TFY", col = "red", breaks = c(0:10), xlab="Number of DB sites (bins)", ylab ="Gene frequency", ylim = c(0,50))
write.table(as.data.frame(gctfy2), file = "TFYcount", quote = FALSE, sep = "\t", row.names = F)


genecountkmo <- read.table("KMOannotated.TSV", header = T, fill = T, row.names = NULL)
par(mfrow= c(1,1))
gckmo <- table(genecountkmo$geneId)
gckmo
gckmo2 <- as.data.frame(table(genecountkmo$geneId))
gckmo2
hist(gckmo2$Freq, xlab = "KMO", col = "grey", main = "", breaks = 10)
hist(gckmo2$Freq, main = "KMO", col = "blue", breaks = c(0:15), xlab="Number of DB sites (bins)", ylab ="Gene frequency", ylim = c(0,500))
write.table(as.data.frame(gckmo2), file = "KMOcount", quote = FALSE, sep = "\t", row.names = F)

genecountkmy <- read.table("KMYannotated.TSV", header = T, fill = T, row.names = NULL)
par(mfrow= c(1,1))
gckmy <- table(genecountkmy$geneId)
gckmy
gckmy2 <- as.data.frame(table(genecountkmy$geneId))
gckmy2
hist(gckmy2$Freq, xlab = "KMY", col = "grey", main = "", breaks = 10)
hist(gckmy2$Freq, main = "KMY", col = "red", breaks = c(0:80), xlab="Number of DB sites (bins)", ylab ="Gene frequency", ylim = c(0,500))
write.table(as.data.frame(gckmy2), file = "KMYcount", quote = FALSE, sep = "\t", row.names = F)


genecountkfo <- read.table("KFOannotated.TSV", header = T, fill = T, row.names = NULL)
par(mfrow= c(1,1))
gckfo <- table(genecountkfo$geneId)
gckfo
gckfo2 <- as.data.frame(table(genecountkfo$geneId))
gckfo2
hist(gckfo2$Freq, xlab = "KFO", col = "grey", main = "", breaks = 10)
hist(gckfo2$Freq, main = "KFO", col = "blue", breaks = c(0:30), xlab="Number of DB sites (bins)", ylab ="Gene frequency", ylim = c(0,500))
write.table(as.data.frame(gckfo2), file = "KFOcount", quote = FALSE, sep = "\t", row.names = F)


genecountkfy <- read.table("KFYannotated.TSV", header = T, fill = T, row.names = NULL)
par(mfrow= c(1,1))
gckfy <- table(genecountkfy$geneId)
gckfy
gckfy2 <- as.data.frame(table(genecountkfy$geneId))
gckfy2
hist(gckfy2$Freq, xlab = "KFY", col = "grey", main = "", breaks = 10)
hist(gckfy2$Freq, main = "KFY", col = "red", breaks = c(0:50), xlab="Number of DB sites (bins)", ylab ="Gene frequency", ylim = c(0,250))
write.table(as.data.frame(gckfy2), file = "KFYcount", quote = FALSE, sep = "\t", row.names = F)





