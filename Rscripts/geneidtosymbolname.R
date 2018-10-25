library(biomaRt)
library(clusterProfiler)
library(org.Hs.eg.db)
library(mygene)



KMOYpromoter <- read.table("promoterKMOY.csv", header = T)
KMOYpromoter
require(clusterProfiler)
#getBM(attributes = c("entrezgene", "external_gene_name", "hgnc_symbol"), filters = "entrezgene", values = TMOgenes, mart = mart)
kashmirPromoter <- bitr(KMOYpromoter$geneId, fromType = "ENTREZID", toType = c("GENENAME", "SYMBOL", "GO" ), OrgDb= "org.Hs.eg.db", drop = T)
kashmirPromoter
write.table(as.data.frame(kashmirPromoter),file = "promogeneinfo2", quote=FALSE, sep="\t", row.names =F )
getGO
#TMO <- read.table("TMOcount", header = T, sep = "\t")
#TMOmart <- bitr(TMO$GeneId, fromType = "ENTREZID", toType=c("GENENAME","SYMBOL", ), OrgDb= "org.Hs.eg.db", drop= FALSE)
sessionInfo()
test <- groupGO(genelist, 'org.Hs.eg.db', keytype = "ENTREZID", ont = "BP", level = 2, readable = F)
genelist <- as.character(KMOYpromoter$geneId)
test
head(summary(test))
write.table(summary(test), file = "test", quote=FALSE, sep="\t", row.names =F )

test2 <- read.table("promoterKMOY.csv", header = T)
test2
res <- queryMany(test2, scopes = 'entrezgene', fields=c('genesymbol', 'go'), species ='human')
res
res[3, 'go.BP'][[3]]
