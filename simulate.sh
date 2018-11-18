#!/bin/bash

s1=$(( $1 - $3 ))
s2=$(( $2 - $3 ))
ovlp=$3

#kg=/lustre/scratch113/teams/zeggini/ReferenceData/1kg/Phase3_May.15/ALL.chr11.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
#kg=/nfs/humgen01/projects/t144_usoc/Coreexome/v1.0/Final_QCed_data/UKHLS_coreex_usgwas_20131119.gencall.smajor_strandupdated_PARupdated_9965samples_525314SNPs_QC
kg=/lustre/scratch114/projects/helic/checkpoints/metacarpa_sim/UKHLS_coreex_usgwas_20131119.gencall.smajor_strandupdated_PARupdated_9965samples_525314SNPs_QC
plink=/software/team144/plink-versions/plink-dev1405/plink

echo "Simulating $s1 and $s2 individuals with $ovlp overlap."

cut -f1 -d' ' $kg.fam | shuf > samples.rnd
head -n $s1 samples.rnd > group1
tail -n +$(( $s1 + 1 )) samples.rnd | head -n $s2 > group2
tail -n +$(( $s1 + $(( $s2 + 1 )) )) samples.rnd | head -n $ovlp > overlap

Rscript --vanilla <(echo "as.matrix(rnorm($(( $(( $s1 + $ovlp )) + $s2 ))))")| awk '{print $2}' |tail -n+2 > pheno

$plink --memory 1000 --bfile $kg --snps-only --keep <(paste <(cat group1 overlap) <(cat group1 overlap)) --pheno <(paste <(cat group1 overlap) <(cat group1 overlap) <(cat <(head -n $s1 pheno) <(tail -n $ovlp pheno)))  --allow-no-sex --out group1 --assoc --make-bed --maf 0.005
sed -r -i 's/^ +//;s/ +/\t/g' group1.qassoc

$plink --memory 1000 --bfile $kg --snps-only --keep <(paste <(cat group2 overlap) <(cat group2 overlap)) --pheno <(paste <(cat group2 overlap) <(cat group2 overlap) <(cat <(head -n $(( $s1 + $s2 )) pheno | tail -n $s2) <(tail -n $ovlp pheno))) --allow-no-sex --out group2 --assoc --make-bed --maf 0.005
sed -r -i 's/^ +//;s/ +/\t/g' group2.qassoc

