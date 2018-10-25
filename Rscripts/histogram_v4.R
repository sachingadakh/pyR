G600 <- read.table("merged.tmp", header = T, fill = T, row.names = NULL, sep = "\t")
head(G600)
#G600[,"H3K9me3_KFO6.Peak.Width"]
names <- colnames(G600)
#count = 0
file=paste("KashmirPeakWidth",".pdf",sep = "")
par(mfrow=c(4,5))
dev.copy(pdf, file)

for (name in names){
  if ((grepl("_K", name)) && (grepl("Peak.Width$", name))) {
    # if (count %% 9 == 0){
    #   dev.off()
    #   file=paste(count,".pdf",sep = "")
    #   par(mfrow=c(3,3))
    #   dev.copy(pdf, file)
    # }
    #count = count + 1
    #print(name)
    hist(G600[,name], xlab = name, col = "grey", main = "", breaks = 10000, xlim = c(0,10000), ylim = c(0,8000))
  }
}
dev.off()

names <- colnames(G600)
#count = 0
file=paste("TamilnaduPeakWidth",".pdf",sep = "")
par(mfrow=c(4,5))
dev.copy(pdf, file)

for (name in names){
  if ((grepl("_T", name)) && (grepl("Peak.Width$", name))) {
    # if (count %% 9 == 0){
    #   dev.off()
    #   file=paste(count,".pdf",sep = "")
    #   par(mfrow=c(3,3))
    #   dev.copy(pdf, file)
    # }
    #count = count + 1
    #print(name)
    hist(G600[,name], xlab = name, col = "grey", main = "", breaks = 10000, xlim = c(0,10000), ylim = c(0,6000))
  }
}
dev.off()

