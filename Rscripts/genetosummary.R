library(rentrez)
library(org.Hs.eg.db)
library(clusterProfiler)

entrez_dbs("gene")
entrez_db_summary("gene")
entrez_db_searchable("gene")
r_search <- entrez_search(db="gene", term = "5629")
r_search
r_search$ids

genetosummary <- read.table("promoterKMOY.csv")
genetosummary
test1 <- entrez_summary(db = "gene", id = genetosummary$V1)
test1
test2 <- extract_from_esummary(test1, c("uid", "nomenclaturename", "summary"))
test2
test3 <- t(test2)
test3
write.table(test3, file = "genetosummary", sep = "\t",  quote = T, row.names = F)


#for (a in genetosummary$V1)
#{test1 <- entrez_summary(db = "gene", id = a)
#id <- a
#sum <- test1$summary
#d=data.frame(id,sum)
#write.table(d, "geneKasummary.csv", sep = "\t", append= TRUE, row.names= FALSE, col.names= FALSE)
#}
#test1
#test1$`284656`$summary
#test1

