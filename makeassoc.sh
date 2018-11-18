#!/bin/bash
heritability=$1
name=betaallequal.$1
/software/bin/Rscript --vanilla /lustre/scratch114/projects/helic/checkpoints/metacarpa_sim/polygeny/makepheno.R $heritability
sleep 5
perl -pe 's/_urn.*?\t/\t/g;s/urn.wtsi./urn:wtsi:/g' phenotype.$name.txt  > phenotype.$name.correctids.txt
/software/team144/plink-versions/beta3v/plink --bfile ../UKHLS_coreex_usgwas_20131119.gencall.smajor_strandupdated_PARupdated_9965samples_525314SNPs_QC --pheno phenotype.$name.correctids.txt --assoc --out $name
awk '{print $NF}' $name.qassoc > $name.pvals
join -1 1 -2 2 <(cut -f2 10k.rep0.traw | sort) <(sort -k2,2 $name.qassoc) > $name.qassoc.onlycausal
join -v 1 -1 2 -2 1 <(sort -k2,2 $name.qassoc) <(cut -f2 10k.rep0.traw | sort) | grep -v NA > $name.qassoc.only_noncausal
cut -d' ' -f1 betaallequal.0.05.qassoc.onlycausal > tagging.$heritability
/software/team144/plink-versions/beta3v/plink --bfile ../UKHLS_coreex_usgwas_20131119.gencall.smajor_strandupdated_PARupdated_9965samples_525314SNPs_QC --show-tags tagging.$heritability --tag-r2 0.2 --tag-kb 500 --out tagging.$heritability
join -v 1 -1 2 -2 1 <(sort -k2,2 $name.qassoc) <(sort tagging.$heritability.tags) | grep -v NA > $name.qassoc.only_noncausal_tags

