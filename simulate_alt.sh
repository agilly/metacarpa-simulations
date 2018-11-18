#!/bin/bash

s1=$(( $1 - $3 ))
s2=$(( $2 - $3 ))
ovlp=$3
maf=$4
effect=$5

#kg=/lustre/scratch113/teams/zeggini/ReferenceData/1kg/Phase3_May.15/ALL.chr11.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
#kg=/nfs/humgen01/projects/t144_usoc/Coreexome/v1.0/Final_QCed_data/UKHLS_coreex_usgwas_20131119.gencall.smajor_strandupdated_PARupdated_9965samples_525314SNPs_QC
kg=/lustre/scratch114/projects/helic/checkpoints/metacarpa_sim/UKHLS_coreex_usgwas_20131119.gencall.smajor_strandupdated_PARupdated_9965samples_525314SNPs_QC
plink=/software/team144/plink-versions/beta3u/plink

echo "Simulating $s1 and $s2 individuals with $ovlp overlap."

cut -f1 -d' ' $kg.fam | shuf > samples.rnd
head -n $s1 samples.rnd > group1
tail -n +$(( $s1 + 1 )) samples.rnd | head -n $s2 > group2
tail -n +$(( $s1 + $(( $s2 + 1 )) )) samples.rnd | head -n $ovlp > overlap

$plink --memory 1000 --bfile $kg --snps-only --keep <(paste <(cat group1 overlap) <(cat group1 overlap)) --allow-no-sex --out group1 --maf 0.001 --make-bed
$plink --memory 1000 --bfile group1 --freq --out group1

$plink --memory 1000 --bfile $kg --snps-only --keep <(paste <(cat group2 overlap) <(cat group2 overlap)) --allow-no-sex --out group2 --maf 0.001 --make-bed
$plink --memory 1000 --bfile group2 --freq --out group2

join -j2 <(sort -k2,2 group1.frq) <(sort -k2,2 group2.frq) | awk -v f=$maf '$5>f*0.95 && $5<f$1.05 && $10>f*0.95 && $10<f$1.05' | head -1 > causal.snp

# Group 1

/software/team144/plink-versions/beta3u/plink --bfile group1 --snp $(cut -f1 -d' ' causal.snp) --recode 01 --output-missing-genotype '.'
count=$(awk 'OFS="\t"{print $1, $1, $8+$7}' plink.ped | cut -f3 |paste -sd+ | bc)

if [ $count -gt $s1 ]
	then
		awk 'OFS="\t"{print $1, $1, 2-($8+$7)}' plink.ped > pheno.g1
	else
		awk 'OFS="\t"{print $1, $1, $8+$7}' plink.ped > pheno.g1
	fi

Rscript --vanilla <(echo 'd=read.table("pheno.g1");d$V3[d$V3!=0]='$effect'*d$V3[d$V3!=0]+rnorm(length(d$V3[d$V3!=0]));d$V3[d$V3==0]=rnorm(length(d$V3[d$V3==0]));options(width=1000);d;') |cut -d' ' -f2- | tail -n+2 | sed -r 's/^\s+//;s/ +/ /g' > pheno.g1f

$plink --memory 1000 --bfile group1 --snps-only --keep <(paste <(cat group1 overlap) <(cat group1 overlap)) --pheno pheno.g1f --allow-no-sex --out group1.pheno --assoc --make-bed

grep -v NA group1.pheno.qassoc > t
mv t group1.pheno.qassoc

rm group1.bed group1.bim group1.fam plink.* 

# Group 2

/software/team144/plink-versions/beta3u/plink --bfile group2 --snp $(cut -f1 -d' ' causal.snp) --recode 01 --output-missing-genotype '.'
count=$(awk 'OFS="\t"{print $1, $1, $8+$7}' plink.ped | cut -f3 |paste -sd+ | bc)

if [ $count -gt $s1 ]
	then
		awk 'OFS="\t"{print $1, $1, 2-($8+$7)}' plink.ped > pheno.g2
	else
		awk 'OFS="\t"{print $1, $1, $8+$7}' plink.ped > pheno.g2
	fi

Rscript --vanilla <(echo 'd=read.table("pheno.g1");d$V3[d$V3!=0]='$effect'*d$V3[d$V3!=0]+rnorm(length(d$V3[d$V3!=0]));d$V3[d$V3==0]=rnorm(length(d$V3[d$V3==0]));options(width=1000);d;') |cut -d' ' -f2- | tail -n+2 | sed -r 's/^\s+//;s/ +/ /g' > pheno.g2f
$plink --memory 1000 --bfile group2 --snps-only --keep <(paste <(cat group1 overlap) <(cat group1 overlap)) --pheno pheno.g2f --allow-no-sex --out group2.pheno --assoc --make-bed
grep -v NA group2.pheno.qassoc > t
mv t group2.pheno.qassoc

rm group2.bed group2.bim group2.fam plink.* 

gzip group1.pheno.bed
gzip group2.pheno.bed
gzip group1.pheno.bim
gzip group2.pheno.bim
rm *.frq
mv group1.pheno.qassoc group1.qassoc
mv group2.pheno.qassoc group2.qassoc


sed -r -i 's/^ +//;s/ +/\t/g' group1.qassoc
sed -r -i 's/^ +//;s/ +/\t/g' group2.qassoc

