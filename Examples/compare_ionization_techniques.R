library(ChEBIplot)

APCI<-read.table('APCI_10k.txt', header=FALSE, sep=' ', nrows=-1)
ESI<-read.table('ESI_10k.txt', header=FALSE, sep=' ', nrows=-1)
EI<-read.table('EI_10k.txt', header=FALSE, sep=' ', nrows=-1)

png('APCI_ESI_EI.png', width=1024, height=1024)
plot3(APCI, ESI, EI, 4, TRUE, ymax=1600, blur=TRUE, tfidf=TRUE)
dev.off()
