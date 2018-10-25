library(ggplot2)
library(reshape2)
library(scales)
library(RColorBrewer)

dbrrepeats <- read.delim("AllDbrRepeats.csv")
head(dbrrepeats)
dbrrepeated<- melt(dbrrepeats)
dbrrepeated
ggplot(data = dbrrepeated, aes(x = variable, y = value, fill = factor(TypeOfRepeats), ordered(TypeOfRepeats)) + geom_bar(stat = "identity") + labs(x="Groups", y="% of repeats") )

dbrrepeats2 <- read.delim("TMO_0.zerocounts", header = T, sep = "\t", stringsAsFactors = FALSE)
dbrrepeats2
ggplot(dbrrepeats2, aes(x=factor(1), y = Frequency, fill=as.factor(paste(NameOfRepeats,Frequency, sep = "-")))) + geom_bar(stat = "identity", width = 1) + ggtitle("TMO:Repeats:1980")+coord_polar(theta = "y")
dbrrepeats2$y = dbrrepeats2$Frequency/2 + c(0, cummax(dbrrepeats2$Frequency) [-length(dbrrepeats2$Frequency)])
ggplot(dbrrepeats2, aes(x=factor(1), y = Frequency, fill=as.factor(paste(NameOfRepeats,Frequency, sep=" - ")))) + geom_bar(stat = "identity", width = 1) + ggtitle("Repeats of TMO:1980") + coord_polar(theta = "y") + xlab("")+ylab("")+theme(legend.position="right", legend.title=element_blank(), plot.title = element_text(lineheight = 3, face = "bold", colour = "black", size = 14)) + geom_text(x = c(1, 1, 1, 1, 1.2, 1.3, 1.4, 1.5), y= Frequency/2 + c(0, cumsum(Frequency)[-length(Frequency)]), label = labels.Frequency)  
                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
ggplot(dbrrepeats2, aes(x=factor(1), y = Frequency, fill=as.factor(paste(NameOfRepeats,Frequency, sep=" - ")))) + geom_bar(stat = "identity", width = 1) + ggtitle("Repeats of TMO:1980") + coord_polar(theta = "y") + xlab("")+ylab("")+theme(legend.position="right", legend.title=element_blank(), plot.title = element_text(lineheight = 3, face = "bold", colour = "black", size = 14)) + geom_text(aes(x = c(1, 1, 1, 1, 1.2, 1.3, 1.4, 1.5), y= Frequency/2 + c(0, cumsum(Frequency)[-length(Frequency)]), label = Frequency), size=3)

ggplot(dbrrepeats2, aes(x=factor(1), y = Frequency, fill=as.factor(paste(NameOfRepeats,Frequency, sep=" - ")))) + geom_bar(stat = "identity", width = 1) + ggtitle("Repeats of TMO:1980") + coord_polar(theta = "y") + xlab("")+ylab("")+theme(legend.position="none", legend.title=element_blank(), plot.title = element_text(lineheight = 3, face = "bold", colour = "black", size = 14)) + geom_text(aes(x = c(1, 1, 1, 1, 1.2, 1.3, 1.4, 1.5, 1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13), y= Frequency/2 + c(0, cumsum(Frequency)[-length(Frequency)]), label = Frequency), size=3)

ggplot(dbrrepeats2, aes(x=factor(1), y = Frequency, fill=as.factor(paste(NameOfRepeats,Frequency, sep=" - ")))) + geom_bar(stat = "identity", width = 1) + ggtitle("Repeats of TMO:1980") + coord_polar(theta = "y") + xlab("")+ylab("")+theme(legend.position="right", legend.title=element_blank(), plot.title = element_text(lineheight = 3, face = "bold", colour = "black", size = 14)) + geom_text(aes(x=factor(1), y=y, label = Frequency), size=3)

