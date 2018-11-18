args=commandArgs(TRUE)
library(data.table)
d=fread(args[1])
ret=0
ret[1]=nrow(d[d$V12<0.05,])/nrow(d)
ret[2]=nrow(d[d$V11<0.05,])/nrow(d)
ret[3]=nrow(d[d$V10<0.05,])/nrow(d)
print(ret)
