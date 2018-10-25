library(org.Hs.eg.db)
library(clusterProfiler)
KFO<- read.table("KFOcount", header=T, sep="\t")
eg.KFO<- bitr(KFO$GeneId, fromType="ENTREZID", toType=c("GENENAME","SYMBOL" ), OrgDb="org.Hs.eg.db", drop= FALSE)
write.table(eg.KFO, "IDsymKFO.txt", sep="\t", quote=F, row.names=F)

KFY<- read.table("KFYcount", header=T, sep="\t")
eg.KFY<- bitr(KFY$GeneId, fromType="ENTREZID", toType=c("GENENAME","SYMBOL" ), OrgDb="org.Hs.eg.db", drop= FALSE)
write.table(eg.KFY, "IDsymKFY.txt", sep="\t", quote=F, row.names=F)

KMO<- read.table("KMOcount", header=T, sep="\t")
eg.KMO<- bitr(KMO$GeneId, fromType="ENTREZID", toType=c("GENENAME","SYMBOL" ), OrgDb="org.Hs.eg.db", drop= FALSE)
write.table(eg.KMO, "IDsymKMO.txt", sep="\t", quote=F, row.names=F)

KMY<- read.table("KMYcount", header=T, sep="\t")
eg.KMY<- bitr(KMY$GeneId, fromType="ENTREZID", toType=c("GENENAME","SYMBOL" ), OrgDb="org.Hs.eg.db", drop= FALSE)
write.table(eg.KMY, "IDsymKMY.txt", sep="\t", quote=F, row.names=F)

TMO<- read.table("TMOcount", header=T, sep="\t")
eg.TMO<- bitr(TMO$GeneId, fromType="ENTREZID", toType=c("GENENAME","SYMBOL" ), OrgDb="org.Hs.eg.db", drop= FALSE)
write.table(eg.TMO, "IDsymTMO.txt", sep="\t", quote=F, row.names=F)

TMY<- read.table("TMYcount", header=T, sep="\t")
eg.TMY<- bitr(TMY$GeneId, fromType="ENTREZID", toType=c("GENENAME","SYMBOL" ), OrgDb="org.Hs.eg.db", drop= FALSE)
write.table(eg.TMY, "IDsymTMY.txt", sep="\t", quote=F, row.names=F)

TFO<- read.table("TFOcount", header=T, sep="\t")
eg.TFO<- bitr(TFO$GeneId, fromType="ENTREZID", toType=c("GENENAME","SYMBOL" ), OrgDb="org.Hs.eg.db", drop= FALSE)
write.table(eg.TFO, "IDsymTFO.txt", sep="\t", quote=F, row.names=F)

TFY<- read.table("TFYcount", header=T, sep="\t")
eg.TFY<- bitr(TFY$GeneId, fromType="ENTREZID", toType=c("GENENAME","SYMBOL" ), OrgDb="org.Hs.eg.db", drop= FALSE)
write.table(eg.TFY, "IDsymTFY.txt", sep="\t", quote=F, row.names=F)
