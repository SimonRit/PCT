#!/bin/sh -x

OPTIONS='--lowmem
    --dimension 700,1,700
    --spacing 1
    --pad 0.
    --hann 0.
    -r proj.*.mhd
    --verbose'

/Users/arbor/Software/rtk/RTK_install/bin/rtksimulatedgeometry \
    -n 360 \
    -f 0 \
    -a -360 \
    -o geometry.rtk \
    --sdd 2400 \
    --sid 2000

/Users/arbor/Workspace/pCT/pct_build/pctfdk \
    ${OPTIONS} \
    -g geometry.rtk \
    -p . \
    -o fdk.mha

$(dirname $0)/make_in_out_proj.sh

/Users/arbor/Software/rtk/RTK_install/bin/rtkfdk \
    ${OPTIONS} \
    -g geometry.rtk \
    -p in  \
    -o fdk_in.mha

/Users/arbor/Software/rtk/RTK_install/bin/rtkfdk \
    ${OPTIONS} \
    -g geometry.rtk \
    -p out \
    -o fdk_out.mha


