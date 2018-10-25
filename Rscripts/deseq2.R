library(DESeq2)
library(GenomicRanges)
library(org.Dm.eg.db)
library(GenomicAlignments)
library(GenomicFeatures)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)

#--------------------------------------------DEseq2_approach_B----------------------------

deseq2 <- read.table("DMbamtowindowcount", header = T, sep = "\t")
rownames(deseq2) <- paste(deseq2$Chr, deseq2$Start, deseq2$Stop, sep = "-")
head(deseq2)
colData <- DataFrame(condition=factor(c("ctrl1", "ctrl1", "treat1", "treat1")))
colData
countData <- as.matrix(deseq2[4:7])
head(countData)
dds <- DESeqDataSetFromMatrix(countData, colData, formula(~ condition))
dds <- DESeq(dds)
res <- results(dds)
res
write.table(res, file = "DEseq2_NRNMR02", sep = "\t", quote = F )

colData2 <- DataFrame(condition=factor(c("ctrl1", "ctrl1", "treat2", "treat2")))
colData2
countData2 <- as.matrix(deseq2[,c(4:5,8:9)])
head(countData2)
dds2 <- DESeqDataSetFromMatrix(countData2, colData2, formula(~ condition))
dds2 <- DESeq(dds2)
res2 <- results(dds2)
res2
write.table(res2, file = "DEseq2_NRNMR1416", sep = "\t", quote = F )

colData3 <- DataFrame(condition=factor(c("ctrl1", "ctrl1", "treat3", "treat3")))
colData3
countData3 <- as.matrix(deseq2[,c(4:5,10:11)])
head(countData3)
dds3 <- DESeqDataSetFromMatrix(countData3, colData3, formula(~ condition))
dds3 <- DESeq(dds3)
res3 <- results(dds3)
res3
write.table(res3, file = "DEseq2_NRNMR016", sep = "\t", quote = F )

 
#-------------------------------------------------------summaryData-----------------------------------


resOrdered <- res[order(res$padj),]
summary(resOrdered$padj)
res2Ordered <- res2[order(res2$padj),]
head(res2Ordered)
summary(res2Ordered$padj)
res3Ordered <- res3[order(res3$padj),]`   `
head(res3Ordered)
summary(res3Ordered$padj)



write.table(resOrdered, file= "padjOfDEseq2NRNM02", sep = "\t", quote = F)
write.table(res2Ordered, file= "padjOfDEseq2NRNM1416", sep = "\t", quote = F)
write.table(res3Ordered, file= "padjOfDEseq2NRNM016", sep = "\t", quote = F)


write.table(resOrdered[1:1000,], file= "top1kpadjOfDEseq2NRNM02", sep = "\t", quote = F)
write.table(res2Ordered[1:1000,], file= "top1kpadjOfDEseq2NRNM1416", sep = "\t", quote = F)
write.table(res3Ordered[1:1000,], file= "top1kpadjOfDEseq2NRNM016", sep = "\t", quote = F)



#--------------------------------------------directCount_approach C----------------------


data1 <- read.table("test016",sep = "\t", header = T)
head(data1)
data1ordere <- data1[rev(order(data1$AvgCount)),]
head(data1ordere)
summary(data1ordere$AvgCount)
write.table(data1ordere[1:1000,], file = "top1kcountNR016", sep = "\t", quote = F, row.names = F)

datanmr02 <- read.table("testnmr02", sep = "\t", header = T)
head(datanmr02)
datanmr02order <- datanmr02[rev(order(datanmr02$AvgCount)),]
head(datanmr02order)
summary(datanmr02order$AvgCount)
write.table(datanmr02order[1:1000,], file = "top1kcountNMR02", sep = "\t", quote = F, row.names = F)

datanmr1416 <- read.table("testnmr1416", sep = "\t", header = T)
head(datanmr1416)
datanmr1416order <- datanmr1416[rev(order(datanmr1416$AvgCount)),]
head(datanmr1416order)
summary(datanmr1416order$AvgCount)
write.table(datanmr1416order[1:1000,], file = "top1kcountNMR1416", sep = "\t", quote = F, row.names = F)

datanmr016 <- read.table("testnmr016", sep = "\t", header = T)
head(datanmr016)
datanmr016order <- datanmr016[rev(order(datanmr016$AvgCount)),]
head(datanmr016order)
summary(datanmr016order$AvgCount)
write.table(datanmr016order[1:1000,], file = "top1kcountNMR016", sep = "\t", quote = F, row.names = F)









