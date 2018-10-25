library(biomaRt)
library(org.Dm.eg.db)

ensembl = useMart("ensembl")
listDatasets(ensembl)
ensembl = useDataset("dmelanogaster_gene_ensembl", mart = ensembl)
listFilters(ensembl)

SigDENRNM02 <- read.table("fbtrsigdiffexpNR16NM2")
head(SigDENRNM02)
SigDENRNM02bcc <- getBM(attributes = c("ensembl_gene_id","ensembl_transcript_id","gene_biotype"), filters = "flybase_transcript_id", values = SigDENRNM02, mart = ensembl)
head(SigDENRNM02bcc)
table(SigDENRNM02bcc$gene_biotype)
write.table(as.data.frame(SigDENRNM02bcc), file= "SigDENRNM02bcc_genebiotype", quote = F, sep = "\t", row.names = F)


SigDENRNM1416 <- read.table("fbtrsigdiffexpNR16NM1416")
head(SigDENRNM1416)
SigDENRNM1416bcc <- getBM(attributes = c("ensembl_gene_id","ensembl_transcript_id","gene_biotype"), filters = "flybase_transcript_id", values = SigDENRNM1416, mart = ensembl)
head(SigDENRNM1416bcc)
table(SigDENRNM1416bcc$gene_biotype)
write.table(as.data.frame(SigDENRNM1416bcc), file= "SigDENRNM1416bcc_genebiotype", quote = F, sep = "\t", row.names = F)


SigDENRNM16 <- read.table("fbtrsigdiffexpNR16NM16")
head(SigDENRNM16)
SigDENRNM16bcc <- getBM(attributes = c("ensembl_gene_id","ensembl_transcript_id","gene_biotype"), filters = "flybase_transcript_id", values = SigDENRNM16, mart = ensembl)
head(SigDENRNM16bcc)
table(SigDENRNM16bcc$gene_biotype)
write.table(as.data.frame(SigDENRNM16bcc), file= "SigDENRNM16bcc_genebiotype", quote = F, sep = "\t", row.names = F)
