args=commandArgs(TRUE)
library(data.table);library(gap)
#file="10k.rep0.traw.0.05.normalised_uniform_standard"
file=args[1]
d=fread(paste(file, ".qassoc.onlycausal", sep=""))
e=fread(paste(file, ".qassoc.only_noncausal", sep=""))
f=fread(paste(file, ".qassoc.only_noncausal_tags", sep=""))
png(paste(file,"graphs.png", sep=""), width=1500, height=1000)
qqunif(as.numeric(e$V9), pch=20)
s=qqunif(as.numeric(d$V9), plot=F)
points(s$x, s$y, pch=20)
s=qqunif(as.numeric(f$V9), plot=F)
points(s$x, s$y, col="forestgreen", pch=20)
dev.off()
