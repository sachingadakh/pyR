library(biomaRt)
library(org.Dm.eg.db)

ensembl = useMart("ensembl")
listDatasets(ensembl)
ensembl = useDataset("dmelanogaster_gene_ensembl", mart = ensembl)
listFilters(ensembl)

NRNM02 <- read.table("fbtrfromtop1knrnm02")
head(NRNM02)
NRNM02bcc <- getBM(attributes = c("ensembl_gene_id","ensembl_transcript_id","gene_biotype"), filters = "flybase_transcript_id", values = NRNM02, mart = ensembl)
head(NRNM02bcc)
table(NRNM02bcc$gene_biotype)
write.table(as.data.frame(NRNM02bcc), file= "NRNM02bcc_genebiotype", quote = F, sep = "\t", row.names = F)

NR016 <- read.table("fbtrtop1kcountnr016")
head(NR016)
NR016bcc <- getBM(attributes = c("ensembl_gene_id","ensembl_transcript_id","gene_biotype"), filters = "flybase_transcript_id", values = NR016, mart = ensembl)
head(NR016bcc)
table(NR016bcc$gene_biotype)
write.table(as.data.frame(NRNM02bcc), file= "NR016bcc_genebiotype", quote = F, sep = "\t", row.names = F)

NMR02 <- read.table("fbtrtop1kcountnmr02")
head(NMR02)
NMR02bcc <- getBM(attributes = c("ensembl_gene_id","ensembl_transcript_id","gene_biotype"), filters = "flybase_transcript_id", values = NMR02, mart = ensembl)
head(NMR02bcc)
table(NMR02bcc$gene_biotype)
write.table(as.data.frame(NMR02bcc), file= "NMR02bcc_genebiotype", quote = F, sep = "\t", row.names = F)

NMR1416 <- read.table("fbtrtop1kcountnmr1416")
head(NMR1416)
NMR1416bcc <- getBM(attributes = c("ensembl_gene_id","ensembl_transcript_id","gene_biotype"), filters = "flybase_transcript_id", values = NMR1416, mart = ensembl)
head(NMR1416bcc)
table(NMR1416bcc$gene_biotype)
write.table(as.data.frame(NMR1416bcc), file= "NMR1416bcc_genebiotype", quote = F, sep = "\t", row.names = F)

NMR016 <- read.table("fbtrtop1kcountnmr016")
head(NMR016)
NMR016bcc <- getBM(attributes = c("ensembl_gene_id","ensembl_transcript_id","gene_biotype"), filters = "flybase_transcript_id", values = NMR016, mart = ensembl)
head(NMR016bcc)
table(NMR016bcc$gene_biotype)
write.table(as.data.frame(NMR016bcc), file= "NMR016bcc_genebiotype", quote = F, sep = "\t", row.names = F)


NR016_RFLTRfbtr <- read.table("uniqFBTRperccNR016")
head(NR016_RFLTRfbtr)
NR016_RFLTRfbtrBCC <- getBM(attributes = c("ensembl_gene_id","ensembl_transcript_id","gene_biotype"), filters = "flybase_transcript_id", values = NR016_RFLTRfbtr, mart = ensembl)
head(NR016_RFLTRfbtrBCC)
table(NR016_RFLTRfbtrBCC$gene_biotype)
cd 


