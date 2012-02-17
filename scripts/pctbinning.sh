#!/bin/sh -x

for NUM in $(seq -w 0 359)
do
    if test -e pairs${NUM}.mhd
    then
#        pctpaircuts \
#          -i pairs${NUM}.mhd \
#          -o pairs${NUM}c.mha \
#          --spacing 2,1.6 \
#          --dimension 190,26 \
#          --source -1000 \
#          --anglecut 3. \
#          --energycut 3.
#        clitkImageInfo pairs${NUM}*.mh?
        pctbinning \
          -i pairs${NUM}.mhd \
          -o proj${NUM}.mhd \
          --spacing=4,1,4 \
          --dimension=95,80,95 \
          --source -1000 \
          --quadricIn="1,0,1,0,0,0,0,0,0,-10000" \
          --count count.mha \
          --anglecut 0.
        clitkImageStatistics -i count.mha -v
#        rm pairs${NUM}c.mha
    fi
done

