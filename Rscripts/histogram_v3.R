G600 <- read.table("merged.tmp", header = T, fill = T, row.names = NULL)
head(G600)

names <- colnames(G600)

file=paste("Kashmir_PeakWidth",".png",sep = "")
par(mfrow=c(1,1))
dev.copy(png, file)
K_PeakWidth <- integer()
for (name in names) {
  if (grepl("Peak.Width$", name )) {
    if (grepl("_K", name )) {
      K_PeakWidth=append(K_PeakWidth,G600[,name])
      # print(name)
      # print(G600[,name])
      
    }
  }
}
hist(K_PeakWidth, xlab = "Kashmir samples Peak width", col = "grey", main = "", breaks = 10000, xlim = c(0,10000))  

dev.off()

file=paste("Tamil_PeakWidth",".png",sep = "")
par(mfrow=c(1,1))
dev.copy(png, file)
T_PeakWidth <- integer()
for (name in names) {
  if (grepl("Peak.Width$", name )) {
    if (grepl("_T", name )) {
      T_PeakWidth=append(T_PeakWidth,G600[,name])
      # print(name)
      # print(G600[,name])
      
    }
  }
}
hist(T_PeakWidth, xlab = "Tamilnadu samples Peak width", col = "grey", main = "", breaks = 10000, xlim = c(0,10000))  

dev.off()
