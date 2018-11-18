#!/bin/bash
gunzip $1/group1.qassoc.gz
gunzip $1/group2.qassoc.gz
/software/bin/Rscript --vanilla getPearsonSingle.R $1
gzip $1/group1.qassoc $1/group2.qassoc
