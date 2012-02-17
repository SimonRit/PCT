#!/bin/sh -x

OPTIONS='-v --lowmem
    --dimension 1024,10,1024
    --spacing 0.25
    --pad 0.0
    -r proj.*.mhd
    --hann 1.'

if ! test -e in
then
    $(dirname $0)/make_in_out_proj.sh
fi

rtksimulatedgeometry \
    -n 360 \
    -f 0 \
    -a -360 \
    -o geometry.rtk \
    --sdd 1600

pctfdk \
    ${OPTIONS} \
    -g geometry.rtk \
    -p . \
    -o fdk.mha

rtkfdk \
    ${OPTIONS} \
    -g geometry.rtk \
    -p in  \
    -o fdk_in.mha

rtkfdk \
    ${OPTIONS} \
    -g geometry.rtk \
    -p out \
    -o fdk_out.mha


rtksimulatedgeometry \
    -n 200 \
    -a -200 \
    -f -320 \
    -o geometry_short.rtk \
    --sdd 1600
rm -fr short
mkdir -p short/in short/out
for i in $(seq -w 320 359) $(seq -w 0 159)
do    
    j=$(echo $i | sed "s/^0*//g")
    if test "$j" = ""
    then
        j=0
    fi
    if test $j -lt 200
    then
        j=$(($j+600))
    fi
    cp proj$i.mhd short/proj$j.mhd
    sed -i "s/ElementDataFile = /ElementDataFile = \.\.\//g" short/proj$j.mhd
    cp in/proj$i.mhd short/in/proj$j.mhd
    sed -i "s/ElementDataFile = /ElementDataFile = \.\.\//g" short/in/proj$j.mhd
    cp out/proj$i.mhd short/out/proj$j.mhd
    sed -i "s/ElementDataFile = /ElementDataFile = \.\.\//g" short/out/proj$j.mhd
done

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

