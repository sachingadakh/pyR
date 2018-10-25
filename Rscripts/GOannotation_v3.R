#library(GenomicRanges)
#library(AnnotationDbi)
#library(GenomicFeatures)
library(GOstats)
library(org.Dr.eg.db)
library("KEGG.db")
#library(tools)

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

inputfile = "./for_venn_diagram_v1/FUG_MUG_v1.txt.csv"

outputfile = "./for_venn_diagram_v1/Female_Control_Up_genes_Ontology_v1.txt"

annotation=org.Dr.egSYMBOL2EG

allgenes = mappedkeys(annotation)
background <- as.list(annotation[allgenes])

query <- read.table(inputfile, sep = "\t")
query = unlist(query)
query <- as.list(annotation[query])

myGOanno(query = query, background = background, annotation = "org.Dr.eg.db", cutoff = 0.01, prefix = 'Up', outfile=outputfile)
