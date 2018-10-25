file1 <- merge(`1htseq`,`2htseq` ,by.x = "V2", by.y= "V2")
file1

df1 <- read.table("htseq1", header = F, sep = "\t")
df2 <- read.table("htseq2", header = F, sep = "\t")
df1
df2
rownames(df1) <- df1$V1
rownames(df2) <- df2$V1
test_df <- merge(df1, df2, by.x = "V1", by.y = "V1")
head(test_df)
write.table(test_df, file = "ht1n2", quote = FALSE, sep = "\t", col.names = T, row.names = F)

df3 <- read.table("htseq3", header = F, sep = "\t")
df4 <- read.table("htseq4", header = F, sep = "\t")
rownames(df3) <- df3$V1
rownames(df4) <- df4$V1
test_df1 <- merge(df3,df4, by.x = "V1", by.y = "V1")
head(test_df1)
write.table(test_df1, file = "ht3n4", quote = F, sep = "\t", col.names = T, row.names = F)

df5 <- read.table("htseq5", header = F, sep = "\t")
df6 <- read.table("htseq6", header = F, sep = "\t")
rownames(df5) <- df5$V1
rownames(df6) <- df6$V1
test_df2 <- merge(df5, df6, by.x = "V1", by.y = "V1")
head(test_df2)
write.table(test_df2, file = "ht5n6", quote = F, sep = "\t", col.names = T, row.names = F)

df7 <- read.table("htseq7", header = F, sep = "\t") 
df8 <- read.table("htseq8", header = F, sep = "\t")
rownames(df7) <- df7$V1
rownames(df8) <- df8$V1
test_df3 <- merge(df7, df8, by.x = "V1", by.y = "V1")
head(test_df3)
write.table(test_df3, file = "ht7n8", quote = F, sep = "\t", col.names = T, row.names = F)
