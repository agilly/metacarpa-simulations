#!/bin/bash
#sleep $(( $4 - 100 ))
mkdir -p $1-$2-$3/$4
cp $6 $1-$2-$3/$4
cd $1-$2-$3/$4
/software/team144/plink-versions/beta3v/plink --memory 5000 --bfile /lustre/scratch114/projects/helic/checkpoints/metacarpa_sim/UKHLS_coreex_usgwas_20131119.gencall.smajor_strandupdated_PARupdated_9965samples_525314SNPs_QC --out causals --maf 0.005 --recode A-transpose --thin-count 20000
/lustre/scratch114/projects/helic/checkpoints/metacarpa_sim/polygeny/simulate_polygeny.sh $1 $2 $3 $5 causals.traw
rm $1-$2-$2/$4/*.bed
rm $1-$2-$2/$4/*.bim
rm $1-$2-$2/$4/*.fam
/nfs/users/nfs_a/ag15/metacapa/linux_64_static_bin/metacapa2 -I group1.qassoc,$1 -I group2.qassoc,$2 -p 9 -b 5 -s 6 -r 2 -c 1 -q 3 -O meta
rm causals.*
gzip group1*
gzip group2*

