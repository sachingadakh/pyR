library(GenomicRanges)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(BSgenome.Hsapiens.UCSC.hg19)
library(org.Hs.eg.db)
library(ggbio) #this is the plotting hero
library(biovizBase)
library(gtools) #for mixedsort()
library(plyr) #revalue() for factors in Events field

#Cytoband
data(ideoCyto, package = "biovizBase")
biovizBase::isIdeogram(ideoCyto$hg19)

p <- autoplot(ideoCyto$hg19, layout="karyogram", cytoband=TRUE) 

#TMO
data <- read.table("TMO_sites", header = T, sep = "\t")
up<-data[data$log2FC >=1 & data$Event == "Up",]
down<-data[data$log2FC <= -1 & data$Event == "Down",]

revalue(down$Event, c("Down" = "Young")) -> down$Event
revalue(up$Event, c("Up" = "Old")) -> up$Event
up[, "strand"] <- NA
down[, "strand"] <- NA
up<- with(up, GRanges(Chrom, IRanges(Start, End), strand, log2FC, Event, pval))
seqlevels(up)

down<- with(down, GRanges(Chrom, IRanges(Start, End), strand, log2FC, Event, pval))
seqlevels(down)

chr.lengths <- seqlengths(Hsapiens)[1:24]
chr.lengths
chr.names <- names(seqlengths(up))
chr.names
chrNames<-mixedsort(chr.names)
chrNames
chr.lengths.match <- chr.lengths[chrNames]
chr.lengths.match
up<- sortSeqlevels(up)
up<- sort(up)

down<- sortSeqlevels(down)
down<- sort(down)

seqlengths(up) <- chr.lengths.match
seqnames(up)
seqlengths(down) <- chr.lengths.match
seqnames(down)
p+ layout_karyogram(up, aes(colour=Event), geom="rect", ylim=c(11,21)) + layout_karyogram(down, aes(colour=Event), geom="rect", ylim=c(11,21)) + scale_colour_manual(values=alpha(c("red", "blue")))


