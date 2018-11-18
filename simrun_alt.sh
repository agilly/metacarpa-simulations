#!/bin/bash
#sleep $(( $4 - 100 ))
mkdir -p $1-$2-$3/$4
cd $1-$2-$3/$4
/lustre/scratch114/projects/helic/checkpoints/metacarpa_sim/code/simulate_alt.sh $1 $2 $3 $5 $6
#rm $1-$2-$2/$4/*.bed
#rm $1-$2-$2/$4/*.bim
#rm $1-$2-$2/$4/*.fam
/nfs/users/nfs_a/ag15/metacapa/metacapa2 -I group1.qassoc,$1 -I group2.qassoc,$2 -p 9 -b 5 -s 6 -r 2 -c 1 -q 3 -O meta
gzip group1.qassoc
gzip group2.qassoc

