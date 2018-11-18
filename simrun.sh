#!/bin/bash
#sleep $(( $4 - 100 ))
mkdir -p $1-$2-$3/$4
cd $1-$2-$3/$4
/nfs/team144/Arthur-Kallia/metacapa_sim/ukhls/simulate.sh $1 $2 $3
rm $1-$2-$2/$4/*.bed
rm $1-$2-$2/$4/*.bim
rm $1-$2-$2/$4/*.fam
/nfs/users/nfs_a/ag15/metacapa/metacapa2 -I group1.qassoc,$1 -I group2.qassoc,$2 -p 9 -b 5 -s 6 -r 2 -c 1 -q 3 -O meta
bgzip group1.qassoc
bgzip group2.qassoc

