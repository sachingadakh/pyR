library(edgeR)
list.files()
cnd1 <- read.table("ht1n2", header = T, sep = "\t")
head(cnd1)
cnd2 <- read.table("ht3n4", header = T , sep = "\t")
head(cnd2)
cnd3 <- read.table("ht5n6", header = T, sep = "\t")
head(cnd3)
cnd4 <- read.table("ht7n8", header = T, sep = "\t")
head(cnd4)

tempdata <- merge(cnd1, cnd2, by="geneid", all = T)
tempdata
head(tempdata)
tempdata <- merge(tempdata, cnd3, by="geneid", all = T)
head(tempdata)
tempdata <- merge(tempdata, cnd4, by="geneid", all = T)
head(tempdata)

#Set first column as row names of dataframe
genes_counts <- data.frame(tempdata[,-1], row.names = tempdata[,1])
head(genes_counts)
write.table(genes_counts, "genescount_allcnd", sep = "\t", quote = F)

genes_cpm <- cpm(genes_counts)
head(genes_cpm)
write.table(genes_cpm, "genesCPM", sep = "\t", quote = F)
