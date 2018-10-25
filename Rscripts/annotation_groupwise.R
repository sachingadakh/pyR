library(GenomicRanges)
library(ChIPseeker)
library(AnnotationDbi)
library(GenomicFeatures)
library(ChIPpeakAnno)
repeats <- read.table("hg19_repeats_selected.bed", sep = "\t", header = T)
repeat_ranges <- toGRanges(repeats, feature="Repeats")
#repeat_ranges <- addMetadata(repeat_ranges, colnames=c("type","names"))
repeat_TxDb <- makeTxDbFromGRanges(repeat_ranges)
repeat_ranges

test_data <- toGRanges(read.table('test.txt', sep = "\t", header = T))
test_data
overlaps.anno <- annotatePeakInBatch(test_data, AnnotationData = repeat_ranges, output = "nearestLocation")
write.csv(as.data.frame(unname(overlaps.anno)), "anno.csv")
pie1(table(overlaps.anno$insideFeature))

TMO <- toGRanges(read.table('TMO.txt', sep = "\t", header = T))
TMO
overlaps.anno <- annotatePeakInBatch(TMO, AnnotationData = repeat_ranges, output = "nearestLocation")

write.csv(as.data.frame(unname(overlaps.anno)), "annoTMO.csv")
pie1(table(overlaps.anno$insideFeature))

TMY <- toGRanges(read.table('TMY.txt', sep = "\t", header = T))
TMY
overlaps.anno <- annotatePeakInBatch(TMY, AnnotationData = repeat_ranges, output = "nearestLocation")
write.csv(as.data.frame(unname(overlaps.anno)), "annoTMY.csv")
pie1(table(overlaps.anno$insideFeature))

for (a in c(TFO,TFY,KMO,KMY,KFO,KFY))
  do 
      for b in $a.txt 
      do          
            $a <- toGRanges(read.table('$b', sep = "\t", header = T));
      
          
      





