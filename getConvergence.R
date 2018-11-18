args=commandArgs(TRUE)
library(methods)
library(Rmpfr)
getcor=function(x,y){return(tetrachoric(table(x,y))$rho);}
getcora=function(x,y){t=table(x,y); alpha=(mpfr(t[1,1], 30)*mpfr(t[2,2], 30))/(mpfr(t[1,2], 30)*mpfr(t[2,1], 30));alpha=alpha^(pi/4);return(as.numeric((alpha-1)/(alpha+1)))}
library(psych)
library(data.table)
d=fread(paste(args[1], "/group1.qassoc", sep=""))
e=fread(paste(args[1], "/group2.qassoc", sep=""))
m=merge(d,e, by="SNP")
range=1:nrow(m)
ma=nrow(m)
numbers=seq(100, ma, by=5000)
o=matrix(0, nrow = length(numbers), ncol = 0);for(i in 1:10){cat(paste("Iteration ", i, "/10\n"));t=NULL;for(number in numbers){cat(paste("Computing with ",number, "SNPs\r"));s=sample(range, number);u=sign(m$BETA.x[s]);v=sign(m$BETA.y[s]);u[u==-1]=0;v[v==-1]=0;t=c(t, getcora(u,v));};cat("\n");o=cbind(o, t)};
write.table(o, paste(gsub("/", "-", args[1]), ".convergence.approx", sep=""), row.names=F, col.names=F, quote=F)
p=matrix(0, nrow = length(numbers), ncol = 0);for(i in 1:10){cat(paste("Iteration ", i, "/10\n"));t=NULL;for(number in numbers){cat(paste("Computing with ",number, "SNPs\r"));s=sample(range, number);u=sign(m$BETA.x[s]);v=sign(m$BETA.y[s]);u[u==-1]=0;v[v==-1]=0;t=c(t, getcor(u,v));};cat("\n");p=cbind(p, t)};
write.table(p, paste(gsub("/", "-", args[1]), ".convergence.ml", sep=""), row.names=F, col.names=F, quote=F) 
