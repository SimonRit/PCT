#!/bin/sh -x

GROUP=20
for FIRST in $(seq 0 $GROUP 1439)
do
qsub -N pctbinning.${FIRST} \
         -o "$(pwd)" \
         -v FIRST=${FIRST},LAST=$(($FIRST+$GROUP-1)),DIRECTORY=$(pwd) \
         /home/srit/src/pct/pct/scripts/pctbinning.job
done

