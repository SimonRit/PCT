#!/bin/sh -x

OPTIONS='--lowmem
    --dimension 210,1,210
    --spacing 1
    --pad 0.
    --hann 0.
    -r proj.*.mhd
    --verbose'

rtksimulatedgeometry \
    -n 360 \
    -f 0 \
    -a -360 \
    -o geometry.rtk \
    --sdd 1110 \
    --sid 1000

pctfdk \
    ${OPTIONS} \
    -g geometry.rtk \
    -p . \
    -o fdk.mha

$(dirname $0)/make_in_out_proj.sh

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


