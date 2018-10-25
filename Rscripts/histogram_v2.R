G600 <- read.table("merged.tmp", header = T, fill = T, row.names = NULL)
head(G600)

names <- colnames(G600)
count = 0
file=paste(count,".png",sep = "")
par(mfrow=c(3,3))
dev.copy(png, file)

for (name in names) {
  if (count %% 9 == 0) {
    dev.off()
    file=paste(count,".png",sep = "")
    par(mfrow=c(3,3))
    dev.copy(png, file)
    
  }
  count = count + 1
  if (grepl("Peak.Width$", name )) {
    hist(G600[,name], xlab = name, col = "grey", main = "", breaks = 10000, xlim = c(0,10000))  
  }
  else if (grepl("FDR$", name)) {
    hist(-log10(G600[,name]), xlab = name, col = "grey", main = "", breaks = 10000, xlim = c(0,50))
  }
  else if (grepl("Fold.Change$", name)) {
    hist(G600[,name], xlab = name, col = "grey", main = "", breaks = 10000, xlim = c(0,10))
  }
}
dev.off()
