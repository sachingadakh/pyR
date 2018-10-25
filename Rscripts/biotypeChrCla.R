library(biomaRt)
library(org.Dm.eg.db)

ensembl = useMart("ensembl")
listDatasets(ensembl)
ensembl = useDataset("dmelanogaster_gene_ensembl",  mart = ensembl)
listFilters(ensembl)
getallgenes <- getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id", "gene_biotype"), mart = ensembl)
head(getallgenes)
write.table(as.data.frame(getallgenes), file = "AllgenesfromDMensembl", quote = F, sep = "\t", row.names = F)
table(getallgenes$gene_biotype)


NRcnd1.1 <- read.table("NRmerged.fpkmabv10")
head(NRcnd1.1)
NRbccabv10 <- getBM(attributes =c("ensembl_gene_id","ensembl_transcript_id","gene_biotype", "entrezgene", "description"), filters = "flybase_transcript_id", values = NRcnd1.1, mart = ensembl)
head(NRbccabv10)
table(NRbccabv10$gene_biotype)
write.table(as.data.frame(NRbccabv10), file = "NRbccabv10_Genebiotype", quote = FALSE, sep = "\t", row.names = F)
write.table(as.data.frame(NRbccabv10), file = "NR016bccabv10_Genebiotypeinbrief", quote = FALSE, sep = "\t", row.names = F)


NRcnd2.1 <- read.table("NRmerged.fpkmblw3", header = F)
NRbccblw3a <- getBM(attributes =c("ensembl_gene_id","ensembl_transcript_id","gene_biotype","entrezgene", "description"), filters = "flybase_transcript_id", values = NRcnd2.1, mart = ensembl)
table(NRbccblw3a$gene_biotype)
head(NRbccblw3a)
write.table(as.data.frame(NRbccblw3a), file = "NRbccblw3_Genebiotype", quote = FALSE, sep = "\t", row.names = F)
write.table(as.data.frame(NRbccblw3a), file = "NR016bccblw3_Genebiotypeinbrief", quote = FALSE, sep = "\t", row.names = F)

NMR02blw3a <- read.table("NMR02merged.fpkmblw3", header = F)
NMR02bccblw3a <- getBM(attributes =c("ensembl_gene_id","ensembl_transcript_id","gene_biotype","entrezgene", "description"), filters = "flybase_transcript_id", values = NMR02blw3a, mart = ensembl)
head(NMR02blw3a)
table(NMR02bccblw3a$gene_biotype)
write.table(as.data.frame(NMR02bccblw3a), file = "NMR02bccblw3_Genebiotype", quote = FALSE, sep = "\t", row.names = F)
write.table(as.data.frame(NMR02bccblw3a), file = "NMR02bccblw3_Genebiotypeinbrief", quote = FALSE, sep = "\t", row.names = F)


NMR02abv10a <- read.table("NMR02merged.fpkmabv10", header = F)
NMR02bccabv10a <- getBM(attributes =c("ensembl_gene_id","ensembl_transcript_id","gene_biotype","entrezgene", "description"), filters = "flybase_transcript_id", values = NMR02abv10a, mart = ensembl)
table(NMR02bccabv10a$gene_biotype)
head(NMR02bccabv10a)
write.table(as.data.frame(NMR02bccabv10a), file = "NMR02bccabv10_Genebiotype", quote = FALSE, sep = "\t", row.names = F)
write.table(as.data.frame(NMR02bccabv10a), file = "NMR02bccabv10_Genebiotypeinbrief", quote = FALSE, sep = "\t", row.names = F)


NMR1416blw3a <- read.table("NMR1416merged.fpkmblw3", header = F)
NMR1416bccblw3a <- getBM(attributes =c("ensembl_gene_id","ensembl_transcript_id","gene_biotype","entrezgene", "description"), filters = "flybase_transcript_id", values = NMR1416blw3a, mart = ensembl)
head(NMR1416bccblw3a)
table(NMR1416bccblw3a$gene_biotype)
write.table(as.data.frame(NMR1416bccblw3a), file = "NMR1416bccblw3_Genebiotype", quote = FALSE, sep = "\t", row.names = F)
write.table(as.data.frame(NMR1416bccblw3a), file = "NMR1416bccblw3_Genebiotypeinbrief", quote = FALSE, sep = "\t", row.names = F)


NMR1416abv10a <- read.table("NMR1416merged.fpkmabv10", header = F)
NMR1416bccabv10a <- getBM(attributes =c("ensembl_gene_id","ensembl_transcript_id","gene_biotype","entrezgene", "description"), filters = "flybase_transcript_id", values = NMR1416abv10a, mart = ensembl)
table(NMR1416bccabv10a$gene_biotype)
head(NMR1416bccabv10a)
write.table(as.data.frame(NMR1416bccabv10a), file = "NMR1416bccabv10_Genebiotype", quote = FALSE, sep = "\t", row.names = F)
write.table(as.data.frame(NMR1416bccabv10a), file = "NMR1416bccabv10_Genebiotypeinbrief", quote = FALSE, sep = "\t", row.names = F)


NMR016blw3a <- read.table("NMR016merged.fpkmblw3", header = F)
NMR016bccblw3a <- getBM(attributes =c("ensembl_gene_id","ensembl_transcript_id","gene_biotype","entrezgene", "description"), filters = "flybase_transcript_id", values = NMR016blw3a, mart = ensembl)
table(NMR016bccblw3a$gene_biotype)
write.table(as.data.frame(NMR016bccblw3a), file = "NMR016bccblw3_Genebiotype", quote = FALSE, sep = "\t", row.names = F)
write.table(as.data.frame(NMR016bccblw3a), file = "NMR016bccblw3_Genebiotypeinbrief", quote = FALSE, sep = "\t", row.names = F)
head(NMR016bccblw3a)

NMR016abv10a <- read.table("NMR016merged.fpkmabv10", header = F)
NMR016bccabv10a <- getBM(attributes =c("ensembl_gene_id","ensembl_transcript_id","gene_biotype","entrezgene", "description"), filters = "flybase_transcript_id", values = NMR016abv10a, mart = ensembl)
table(NMR016bccabv10a$gene_biotype)
head(NMR016bccabv10a)
write.table(as.data.frame(NMR016bccabv10a), file = "NMR016bccabv10_Genebiotype", quote = FALSE, sep = "\t", row.names = F)
write.table(as.data.frame(NMR016bccabv10a), file = "NMR016bccabv10_Genebiotypeinbrief", quote = FALSE, sep = "\t", row.names = F)

NR3t10 <- read.table("NRmerged.fpkm3to10", header = F)
NR016bcc3t10 <- getBM(attributes =c("ensembl_gene_id","ensembl_transcript_id","gene_biotype","entrezgene", "description"), filters = "flybase_transcript_id", values = NR3t10, mart = ensembl)
head(NR016bcc3t10)
write.table(as.data.frame(NR016bcc3t10), file = "NR016bcc3t10_Genebiotype", quote = FALSE, sep = "\t", row.names = F)
write.table(as.data.frame(NR016bcc3t10), file = "NR016bcc3t10_Genebiotypeinbrief", quote = FALSE, sep = "\t", row.names = F)
table(NR016bcc3t10$gene_biotype)

NMR02cnd3t10 <- read.table("NMR02merged.fpkm3to10", header = F)
NMR02bcc3t10 <- getBM(attributes =c("ensembl_gene_id","ensembl_transcript_id","gene_biotype","entrezgene", "description"), filters = "flybase_transcript_id", values = NMR02cnd3t10, mart = ensembl)
table(NMR02bcc3t10$gene_biotype)
head(NMR02bcc3t10)
write.table(as.data.frame(NMR02bcc3t10), file = "NMR02bcc3t10_Genebiotype", quote = FALSE, sep = "\t", row.names = F)
write.table(as.data.frame(NMR02bcc3t10), file = "NMR02bcc3t10_Genebiotypeinbrief", quote = FALSE, sep = "\t", row.names = F)

NMR1416cnd3t10 <- read.table("NMR1416merged.fpkm3to10", header = F)
NMR1416bcc3t10 <- getBM(attributes =c("ensembl_gene_id","ensembl_transcript_id","gene_biotype","entrezgene", "description"), filters = "flybase_transcript_id", values = NMR1416cnd3t10, mart = ensembl)
head(NMR1416bcc3t10)
table(NMR1416bcc3t10$gene_biotype)
write.table(as.data.frame(NMR1416bcc3t10), file = "NMR1416bcc3t10_Genebiotype", quote = FALSE, sep = "\t", row.names = F)
write.table(as.data.frame(NMR1416bcc3t10), file = "NMR1416bcc3t10_Genebiotypeinbrief", quote = FALSE, sep = "\t", row.names = F)

NMR016cnd3t10 <- read.table("NMR016merged.fpkm3to10", header = F)
NMR016bcc3t10 <- getBM(attributes =c("ensembl_gene_id","ensembl_transcript_id","gene_biotype","entrezgene", "description"), filters = "flybase_transcript_id", values = NMR016cnd3t10, mart = ensembl)
head(NMR016bcc3t10)
table(NMR016bcc3t10$gene_biotype)
write.table(as.data.frame(NMR016bcc3t10), file = "NMR016bcc3to10_Genebiotype", quote = FALSE, sep = "\t", row.names = F)
write.table(as.data.frame(NMR016bcc3t10), file = "NMR016bcc3to10_Genebiotypeinbrief", quote = FALSE, sep = "\t", row.names = F)


NRtrsize <- read.table("nr.trsize", header = F)
head(NRtrsize)
mergingNRsizenbiotype <- merge.data.frame(NRbccabv10, NRtrsize,  by.x = "ensembl_transcript_id", by.y = "V1", all.x = T)
head(mergingNRsizenbiotype)
write.table(as.data.frame(mergingNRsizenbiotype), file = "NR016_withtranscriptsize", quote = FALSE, sep = "\t", row.names = F)

NMR02trsize <- read.table("nmr02.trsize", header = F)
head(NMR02trsize)
mergingNMR02sizenbiotype <- merge.data.frame(NMR02bccabv10a, NMR02trsize,  by.x = "ensembl_transcript_id", by.y = "V1", all.x = T)
head(mergingNMR02sizenbiotype)
write.table(as.data.frame(mergingNMR02sizenbiotype), file = "NMR02_withtranscriptsize", quote = FALSE, sep = "\t", row.names = F)


NMR1416trsize <- read.table("nmr1416.trsize", header = F)
head(NMR1416trsize)
mergingNMR1416sizenbiotype <- merge.data.frame(NMR1416bccabv10a, NMR1416trsize,  by.x = "ensembl_transcript_id", by.y = "V1", all.x = T)
head(mergingNMR1416sizenbiotype)
write.table(as.data.frame(mergingNMR1416sizenbiotype), file = "NMR1416_withtranscriptsize", quote = FALSE, sep = "\t", row.names = F)

NMR016trsize <- read.table("nmr016.trsize", header = F)
head(NMR016trsize)
mergingNMR016sizenbiotype <- merge.data.frame(NMR016bccabv10a, NMR016trsize,  by.x = "ensembl_transcript_id", by.y = "V1", all.x = T)
head(mergingNMR016sizenbiotype)
write.table(as.data.frame(mergingNMR016sizenbiotype), file = "NMR016_withtranscriptsize", quote = FALSE, sep = "\t", row.names = F)





