#!/bin/sh -x

OPTIONS='--lowmem
    --dimension 1000,1,1000
    --spacing 0.25
    --pad 0.
    --hann 0.
    -r proj.*.mhd
    --verbose'

rtksimulatedgeometry \
    -n 720 \
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


