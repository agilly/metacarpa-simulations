
args=commandArgs(TRUE)
heritability=args[1]
library(data.table)
d=fread("10k.rep0.traw")
d=data.frame(d)
dprime=d[,7:(ncol(d)-2)]
colnames(dprime)=""
af=rowSums(dprime, na.rm=T)/(2*ncol(dprime))
heritability=0.05
betas=rnorm(nrow(dprime), sd=heritability/nrow(dprime))
gi=apply(dprime, 2, function(x) sum(x*betas, na.rm=T))
sgi=var(gi, na.rm=T)
erterm=rnorm(ncol(dprime), mean=0, sd=sgi*(1-heritability)/heritability)
write.table(cbind(colnames(d)[7:(ncol(d)-2)],colnames(d)[7:(ncol(d)-2)],erterm+gi), paste("phenotype.betaallequal.", args[1],".txt",sep=""), row.names=F, col.names=F, quote=F, sep="\t")
