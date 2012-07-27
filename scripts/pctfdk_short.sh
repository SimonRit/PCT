#!/bin/sh -x

OPTIONS='--lowmem
    --dimension 2100,1,2100
    --spacing 0.1
    --pad 0.
    --hann 0.
    -r proj.*.mhd
    --verbose'

rm -fr short
mkdir -p short/in short/out
for i in $(seq -w 640 719) $(seq -w 0 319)
do    
    j=$(echo $i | sed "s/^0*//g")
    if test "$j" = ""
    then
        j=0
    fi
    if test $j -lt 400
    then
        j=$(($j+600))
    elif test $j -gt 600
    then
        j=$(($j-500))
    fi
    cp proj0$i.mhd short/proj0$j.mhd
    sed -i "s/ElementDataFile = /ElementDataFile = \.\.\//g" short/proj0$j.mhd
    cp in/proj0$i.mhd short/in/proj0$j.mhd
    sed -i "s/ElementDataFile = /ElementDataFile = \.\.\//g" short/in/proj0$j.mhd
    cp out/proj0$i.mhd short/out/proj0$j.mhd
    sed -i "s/ElementDataFile = /ElementDataFile = \.\.\//g" short/out/proj0$j.mhd
done

rtksimulatedgeometry \
    -n 400 \
    -a -200 \
    -f -320 \
    -o geometry_short.rtk \
    --sdd 1110 \
    --sid 1000

pctfdk \
    ${OPTIONS} \
    -g geometry_short.rtk \
    -p short \
    -o fdk_short.mha

rtkfdk \
    ${OPTIONS} \
    -g geometry_short.rtk \
    -p short/in  \
    -o fdk_short_in.mha

rtkfdk \
    ${OPTIONS} \
    -g geometry_short.rtk \
    -p short/out \
    -o fdk_short_out.mha

