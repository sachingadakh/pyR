library(GOstats)
library(org.Dm.eg.db)
library(KEGG.db)

myGOanno <- function (query, background, annotation, cutoff, prefix='GO', outfile) {
  directions = c('over', 'under')
  ontologies = c('BP', 'MF', 'CC')
  for(d in directions) {
    for(o in ontologies) {
      GOparams <- new("GOHyperGParams", geneIds=query, universeGeneIds=background, annotation=annotation, ontology=o, pvalueCutoff=cutoff, conditional=F, testDirection=d)
      cat(sprintf("Performing test for %s in %s direction\n", o, d), file = stderr())
      test <- hyperGTest(GOparams)
      write.table(summary(test), file = outputfile, quote = 7, append = TRUE, row.names = F, sep="\t")
    }
    cat(sprintf("Performing test for KEGG in %s direction\n", d), file = stderr())
    params <- new("KEGGHyperGParams", geneIds=query, universeGeneIds=background, annotation=annotation, pvalueCutoff=cutoff, testDirection=d)
    hgOver <- hyperGTest(params)
    write.table(summary(hgOver), file = outputfile, quote = 7, append = TRUE, row.names = F, sep="\t")
  }
}

inputfile = "idnname"
inputfile
outputfile = "idsnr_GO"
annotation =org.Dm.egSYMBOL2EG
allgenes = mappedkeys(annotation)
background <- as.list(annotation[allgenes])

query <- read.table(inputfile, sep = "\t")
query = unlist(query)
query <- as.list(annotation[query])

myGOanno(query = query, background = background, annotation = "org.Dm.eg.db", cutoff = 0.01, prefix = 'Up', outfile = outputfile)
