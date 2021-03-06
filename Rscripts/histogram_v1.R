G450 <- read.table("W150G450pwfcfdr_merged.txt", header = T, fill = T, row.names = NULL)
head(G450)
par(mfrow=c(3,3))
hist(G450$PeakWidth_G450highest, xlab = "Peak Width W150 G450 FDRe-2 HPN", col = "grey", main = "", breaks=10000, xlim = c(0, 10000))
hist(-log10(G450$FDR_G450highest), xlab = "FDR W150 G450 FDRe-2 HPN", col = "grey", main = "", breaks = 10000, xlim = c(0, 50))
hist(G450$FoldChange_G450highest, xlab = "Fold Change W150 G450 FDRe-2 HPN", col = "grey", main = "", breaks = 10000, xlim = c(0, 10))
hist(G450$PeakWidth_G450median, xlab = "Peak Width 150 G450 FDRe-2 MPN", col = "grey", main = "", breaks =10000, xlim = c(0, 10000))
hist(-log10(G450$FDR_G450median), xlab = "FDR W150 G450 FDRe-2 MPN", col = "grey", main = "",breaks = 10000, xlim = c(0, 50))
hist(G450$FoldChange_G450median, xlab = "Fold Change W150 G450 FDRe-2 MPN", col = "grey", main = "", breaks = 10000, xlim = c(0, 10))
hist(G450$PeakWidth_G450lowest, xlab = "Peak Width W150 G450 FDRe-2 LPN", col = "grey", main = "", breaks=10000, xlim = c(0, 10000))
hist(-log10(G450$FDR_G450lowest), xlab = "FDR W150 G450 FDRe-2 LPN", col = "grey", main = "",breaks = 10000, xlim = c(0, 50))
hist(G450$FoldChange_G450lowest, xlab = "Fold Change W150 G450 FDRe-2 LPN", col = "grey", main = "", breaks = 10000, xlim = c(0, 10))

G400 <- read.table("W200G400pwfcfdr_merged.txt", header = T, fill = T, row.names = NULL)
head(G400)
par(mfrow=c(3,3))
hist(G400$PeakWidth_G400highest, xlab = "Peak Width W200 G400 FDRe-2 HPN", col = "grey", main = "", breaks=10000, xlim = c(0, 10000))
hist(-log10(G400$FDR_G400highest), xlab = "FDR W200 G400 FDRe-2 HPN", col = "grey", main = "", breaks = 10000, xlim = c(0, 50))
hist(G400$FoldChange_G400highest, xlab = "Fold Change W200 G400 FDRe-2 HPN", col = "grey", main = "", breaks = 10000, xlim = c(0, 10))
hist(G400$PeakWidth_G400median, xlab = "Peak Width W200 G400 FDRe-2 MPN", col = "grey", main = "", breaks=10000, xlim = c(0, 10000))
hist(-log10(G400$FDR_G400median), xlab = "FDR W200 G400 FDRe-2 MPN", col = "grey", main = "", breaks = 10000, xlim = c(0, 50))
hist(G400$FoldChange_G400median, xlab = "Fold Change W200 G400 FDRe-2 MPN", col = "grey", main = "", breaks = 10000, xlim = c(0, 10))
hist(G400$PeakWidth_G400lowest, xlab = "Peak Width W200 G400 FDRe-2 LPN", col = "grey", main = "", breaks=10000, xlim = c(0, 10000))
hist(-log10(G400$FDR_G400lowest), xlab = "FDR W200 G400 FDRe-2 LPN", col = "grey", main = "", breaks = 10000, xlim = c(0, 50))
hist(G400$FoldChange_G400lowest, xlab = "Fold Change W200 G400 FDRe-2 LPN", col = "grey", main = "", breaks = 10000, xlim = c(0, 10))

G600 <- read.table("W200G600pwfcfdr_merged.txt", header = T, fill = T, row.names = NULL)
head(G600)
par(mfrow=c(3,3))
hist(G600$PeakWidth_G600highest, xlab = "Peak Width W200 G600 FDRe-2 HPN", col = "grey", main = "", breaks=10000, xlim = c(0, 10000))
hist(-log10(G600$FDR_G600highest), xlab = "FDR W200 G600 FDRe-2 HPN", col = "grey", main = "", breaks = 10000, xlim = c(0, 50))
hist(G600$FoldChange_G600highest, xlab = "Fold Change W200 G600 FDRe-2 HPN", col = "grey", main = "", breaks = 10000, xlim = c(0, 10))
hist(G600$PeakWidth_G600median, xlab = "Peak Width W200 G600 FDRe-2 MPN", col = "grey", main = "", breaks=10000, xlim = c(0, 10000))
hist(-log10(G600$FDR_G600median), xlab = "FDR W200 G600 FDRe-2 MPN", col = "grey", main = "", breaks = 10000, xlim = c(0, 50))
hist(G600$FoldChange_G600median, xlab = "Fold Change W200 G600 FDRe-2 MPN", col = "grey", main = "", breaks = 10000, xlim = c(0, 10))
hist(G600$PeakWidth_G600lowest, xlab = "Peak Width W200 G600 FDRe-2 LPN", col = "grey", main = "", breaks=10000, xlim = c(0, 10000))
hist(-log10(G600$FDR_G600lowest), xlab = "FDR W200 G600 FDRe-2 LPN", col = "grey", main = "", breaks = 10000, xlim = c(0, 50))
hist(G600$FoldChange_G600lowest, xlab = "Fold Change W200 G600 FDRe-2 LPN", col = "grey", main = "", breaks = 10000, xlim = c(0, 10))


