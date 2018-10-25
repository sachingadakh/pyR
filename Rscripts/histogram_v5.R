G600_new <- read.table("W200_G600_peak_widths.txt", header = T, row.names = NULL, sep = "\t", fill = T)
head(G600_new)

names <- colnames(G600_new)
# library(compare)
# comparison <- compare(G600_new[,"KFO6"],G600[,"H3K9me3_KFO6.Peak.Width"],allowAll=TRUE)
# comparison$tM
# a1<-as.data.frame(G600_new[,"KFO6"])
# head(a1)
# nrow(a1)
# a2<-as.data.frame(G600[,"H3K9me3_KFO6.Peak.Width"])
# head(a2)
# nrow(a2)
# difference <- data.frame(lapply(1:ncol(G600_new[,"KFO6"]),function(i)setdiff(G600_new[,"KFO6"][,i],comparison$tM[,i])))
# write.table(a1,"a1.txt")
# write.table(a2,"a2.txt")
# require(sqldf)
# a1NotIna2 <- sqldf('SELECT * FROM a2 EXCEPT SELECT * FROM a1')
# head(a1NotIna2)
#count = 0

file=paste("KashmirPeakWidth","_v1.pdf",sep = "")
par(mfrow=c(1,1))
dev.copy(pdf, file)
K_PeakWidth <- integer()
for (name in names) {
    if (grepl("^K", name )) {
      K_PeakWidth=append(K_PeakWidth,G600_new[,name])
      # print(name)
      # print(G600[,name])
      
    }
}
hist(K_PeakWidth, xlab = "Kashmir samples Peak width", col = "grey", main = "", breaks = 10000, xlim = c(0,10000))


par(mfrow=c(4,5))


for (name in names){
  if (grepl("^K", name)) {
    # if (count %% 9 == 0){
    #   dev.off()
    #   file=paste(count,".pdf",sep = "")
    #   par(mfrow=c(3,3))
    #   dev.copy(pdf, file)
    # }
    #count = count + 1 G600[,name]
    #print(name)
    hist(G600_new[,name], xlab = name, col = "grey", main = "", breaks = 10000, xlim = c(0,10000), ylim = c(0,8000))
  }
}
dev.off()







file=paste("TamilnaduPeakWidth","_v1.pdf",sep = "")
par(mfrow=c(1,1))
dev.copy(pdf, file)
T_PeakWidth <- integer()
for (name in names) {
  if (grepl("^T", name )) {
    T_PeakWidth=append(T_PeakWidth,G600_new[,name])
    # print(name)
    # print(G600[,name])
    
  }
}
hist(T_PeakWidth, xlab = "Tamilnadu samples Peak width", col = "grey", main = "", breaks = 10000, xlim = c(0,10000))


par(mfrow=c(4,5))


for (name in names){
  if (grepl("^T", name)) {
    # if (count %% 9 == 0){
    #   dev.off()
    #   file=paste(count,".pdf",sep = "")
    #   par(mfrow=c(3,3))
    #   dev.copy(pdf, file)
    # }
    #count = count + 1 G600[,name]
    #print(name)
    hist(G600_new[,name], xlab = name, col = "grey", main = "", breaks = 10000, xlim = c(0,10000), ylim = c(0,8000))
  }
}
dev.off()



