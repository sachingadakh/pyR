library(GenomicRanges)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(BSgenome.Hsapiens.UCSC.hg19)
library(org.Hs.eg.db)
library(ggbio) #this is the plotting hero
library(biovizBase)
library(gtools) #for mixedsort()
library(plyr) #revalue() for factors in Events field

#TM
data <- read.table("TMO_sites", header = T, sep = "\t")
up<-data[data$log2FC >=1 & data$Event == "Up",]
down<-data[data$log2FC <= -0.99 & data$Event == "Down",]	
TMO<-rbind(up,down)
revalue(TMO$Event, c("Down" = "Young")) -> TMO$Event
revalue(TMO$Event, c("Up" = "Old")) -> TMO$Event
TMO[, "strand"] <- NA
TMObed <- with(TMO, GRanges(Chrom, IRanges(Start, End), strand, log2FC, Event, pval))
seqlevels(TMObed)
chr.lengths <- seqlengths(Hsapiens)[1:24]
chr.lengths
chr.names <- names(seqlengths(TMObed))
chr.names
chrNames<-mixedsort(chr.names)
chrNames
chr.lengths.match <- chr.lengths[chrNames]
chr.lengths.match
TMObed<- sortSeqlevels(TMObed)
TMObed<- sort(TMObed)
seqlengths(TMObed) <- chr.lengths.match
seqnames(TMObed)
autoplot(TMObed, layout="karyogram", aes(color = Event, fill = Event))

#TF
data <- read.table("TFO_sites", header = T, sep = "\t")
up<-data[data$log2FC >=1 & data$Event == "Up",]
down<-data[data$log2FC <= -0.99 & data$Event == "Down",]
TFO<-rbind(up,down)
revalue(TFO$Event, c("Down" = "Young")) -> TFO$Event
revalue(TFO$Event, c("Up" = "Old")) -> TFO$Event
TFO[, "strand"] <- NA
TFObed <- with(TFO, GRanges(Chrom, IRanges(Start, End), strand, log2FC, Event, pval))
chr.names <- names(seqlengths(TFObed))
chrNames<-mixedsort(chr.names)
chr.lengths.match <- chr.lengths[chrNames]
TFObed<- sortSeqlevels(TFObed)
TFObed<- sort(TFObed)
seqlengths(TFObed) <- chr.lengths.match
autoplot(TFObed, layout="karyogram", aes(color = Event, fill = Event))

#KM
data <- read.table("KMO_sites", header = T, sep = "\t")
up<-data[data$log2FC >=1 & data$Event == "Up",]
down<-data[data$log2FC <= -0.99 & data$Event == "Down",]
KMO<-rbind(up,down)
revalue(KMO$Event, c("Down" = "Young")) -> KMO$Event
revalue(KMO$Event, c("Up" = "Old")) -> KMO$Event
KMO[, "strand"] <- NA
KMObed <- with(KMO, GRanges(Chrom, IRanges(Start, End), strand, log2FC, Event, pval))
seqlevels(KMObed)
chr.lengths <- seqlengths(Hsapiens)[1:24]
chr.names <- names(seqlengths(KMObed))
chrNames<-mixedsort(chr.names)
chr.lengths.match <- chr.lengths[chrNames]
KMObed<- sortSeqlevels(KMObed)
KMObed<- sort(KMObed)
seqlengths(KMObed) <- chr.lengths.match
seqnames(KMObed)
autoplot(KMObed, layout="karyogram", aes(color = Event, fill = Event))

#KF
data <- read.table("KFO_sites", header = T, sep = "\t")
up<-data[data$log2FC >=1 & data$Event == "Up",]
down<-data[data$log2FC <= -0.99 & data$Event == "Down",]
KFO<-rbind(up,down)
revalue(KFO$Event, c("Down" = "Young")) -> KFO$Event
revalue(KFO$Event, c("Up" = "Old")) -> KFO$Event
KFO[, "strand"] <- NA
KFObed <- with(KFO, GRanges(Chrom, IRanges(Start, End), strand, log2FC, Event, pval))
chr.names <- names(seqlengths(KFObed))
chrNames<-mixedsort(chr.names)
chr.lengths.match <- chr.lengths[chrNames]
KFObed<- sortSeqlevels(KFObed)
KFObed<- sort(KFObed)
seqlengths(KFObed) <- chr.lengths.match
autoplot(KFObed, layout="karyogram", aes(color = Event, fill = Event))
