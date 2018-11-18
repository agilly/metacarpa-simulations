args=commandArgs(TRUE)
print(paste("Got file", args[1], " and heritability ", args[2]))
heritability=as.numeric(args[2])
library(data.table)
d=fread(args[1])
d=data.frame(d)
dprime=d[,7:(ncol(d)-2)]
colnames(dprime)=""
af=rowSums(dprime, na.rm=T)/(2*ncol(dprime))
erf.inv <- function(x) qnorm((x + 1)/2)/sqrt(2)
#betas=rnorm(nrow(dprime), sd=heritability/nrow(dprime))
betas=sample(c(-1,1), length(af), replace=T)*(0.5-(af-0.005)/(1-0.005))
betas=sqrt(2)*erf.inv(2*betas)*(heritability/length(betas))

mu_g=2*af
sigma_g=sqrt(2*af*(1-af))
z=apply(dprime, 2, function(x) (x-mu_g)/sigma_g)
gi=apply(z, 2, function(x) sum(x*betas, na.rm=T))
sgi=var(gi, na.rm=T)
erterm=rnorm(ncol(dprime), mean=0, sd=sgi*(1-heritability)/heritability)

#gi=apply(dprime, 2, function(x) sum(x*betas, na.rm=T))
#sgi=var(gi, na.rm=T)
#erterm=rnorm(ncol(dprime), mean=0, sd=sgi*(1-heritability)/heritability)
write.table(cbind(colnames(d)[7:(ncol(d)-2)],colnames(d)[7:(ncol(d)-2)],erterm+gi), paste("phenotype.", args[1], ".", heritability,".normalised_uniform_standard.txt",sep=""), row.names=F, col.names=F, quote=F, sep="\t")
