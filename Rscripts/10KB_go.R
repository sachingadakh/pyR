library(GenomicRanges)
library(AnnotationDbi)
library(GenomicFeatures)
library(GOstats)
library(org.Hs.eg.db)
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
inputfile = "KFO_10kb"
outputfile = "KFO10kb_GO"
annotation =org.Hs.egSYMBOL2EG
allgenes = mappedkeys(annotation)
background <- as.list(annotation[allgenes])

query <- read.table(inputfile, sep = "\t")
query = unlist(query)
query <- as.list(annotation[query])

myGOanno(query = query, background = background, annotation = "org.Hs.eg.db", cutoff = 0.01, prefix = 'Up', outfile = outputfile)


inputfile = "KFY_10kb"
outputfile = "KFY10kb_GO"
annotation =org.Hs.egSYMBOL2EG
allgenes = mappedkeys(annotation)
background <- as.list(annotation[allgenes])

query <- read.table(inputfile, sep = "\t")
query = unlist(query)
query <- as.list(annotation[query])

myGOanno(query = query, background = background, annotation = "org.Hs.eg.db", cutoff = 0.01, prefix = 'Up', outfile = outputfile)


inputfile = "KMO_10kb"
outputfile = "KMO10kb_GO"
annotation =org.Hs.egSYMBOL2EG
allgenes = mappedkeys(annotation)
background <- as.list(annotation[allgenes])

query <- read.table(inputfile, sep = "\t")
query = unlist(query)
query <- as.list(annotation[query])

myGOanno(query = query, background = background, annotation = "org.Hs.eg.db", cutoff = 0.01, prefix = 'Up', outfile = outputfile)

#inputfile = "KMY_10kb"
outputfile = "KMY10kb_GO"
annotation =org.Hs.egSYMBOL2EG
allgenes = mappedkeys(annotation)
background <- as.list(annotation[allgenes])

query <- read.table(inputfile, sep = "\t")
query = unlist(query)
query <- as.list(annotation[query])

myGOanno(query = query, background = background, annotation = "org.Hs.eg.db", cutoff = 0.01, prefix = 'Up', outfile = outputfile)

inputfile = "TMO_10kb"
outputfile = "TMO10kb_GO"
annotation =org.Hs.egSYMBOL2EG
allgenes = mappedkeys(annotation)
background <- as.list(annotation[allgenes])

query <- read.table(inputfile, sep = "\t")
query = unlist(query)
query <- as.list(annotation[query])

myGOanno(query = query, background = background, annotation = "org.Hs.eg.db", cutoff = 0.01, prefix = 'Up', outfile = outputfile)

inputfile = "TMY_10kb"
outputfile = "TMY10kb_GO"
annotation =org.Hs.egSYMBOL2EG
allgenes = mappedkeys(annotation)
background <- as.list(annotation[allgenes])

query <- read.table(inputfile, sep = "\t")
query = unlist(query)
query <- as.list(annotation[query])

myGOanno(query = query, background = background, annotation = "org.Hs.eg.db", cutoff = 0.01, prefix = 'Up', outfile = outputfile)

inputfile = "TFO_10kb"
outputfile = "TFO10kb_GO"
annotation =org.Hs.egSYMBOL2EG
allgenes = mappedkeys(annotation)
background <- as.list(annotation[allgenes])

query <- read.table(inputfile, sep = "\t")
query = unlist(query)
query <- as.list(annotation[query])

myGOanno(query = query, background = background, annotation = "org.Hs.eg.db", cutoff = 0.01, prefix = 'Up', outfile = outputfile)

inputfile = "TFY_10kb"
outputfile = "TFY10kb_GO"
annotation =org.Hs.egSYMBOL2EG
allgenes = mappedkeys(annotation)
background <- as.list(annotation[allgenes])

query <- read.table(inputfile, sep = "\t")
query = unlist(query)
query <- as.list(annotation[query])

myGOanno(query = query, background = background, annotation = "org.Hs.eg.db", cutoff = 0.01, prefix = 'Up', outfile = outputfile)

inputfile = "promoter.TSV"
outputfile = "promoterKMOY_GO"
annotation =org.Hs.egSYMBOL2EG
allgenes = mappedkeys(annotation)
background <- as.list(annotation[allgenes])

query <- read.table(inputfile, sep = "\t")
query = unlist(query)
query <- as.list(annotation[query])

myGOanno(query = query, background = background, annotation = "org.Hs.eg.db", cutoff = 0.01, prefix = 'Up', outfile = outputfile)
