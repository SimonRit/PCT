#!/bin/sh -x

for FILENAME in *gz
do
    tar xvfz $FILENAME
done

for FILENAME in $(ls */PhaseSpaceIn.root)
do
    cd $(dirname $FILENAME)
    ../$(dirname $0)/pctpairprotons.sh
    rm PhaseSpace*root
    rm *txt
    cd ..
done

for FILENAME in $(ls */pairs000.mhd)
do
    if ! test -e merged
    then
        mkdir merged
        cp $(dirname $FILENAME)/pairs???.mhd merged
    else
        for i in $(dirname $FILENAME)/pairs???.mhd
        do
            CURRENTSIZE=$(grep DimSize merged/$(basename $i) | sed 's/.*DimSize = \(.*\) \(.*\)/\2/g' )
            NEWSIZE=$(grep DimSize $i | sed 's/.*DimSize = \(.*\) \(.*\)/\2/g' )
            SUMSIZE=$(($CURRENTSIZE+$NEWSIZE))
            sed -i "s/.*DimSize = \(.*\) \(.*\)/DimSize = \1 $(($CURRENTSIZE+$NEWSIZE))/g" merged/$(basename $i)
        done
    fi
done

cd merged
for i in pairs???.mhd
do
    cat ../*/$(basename $i mhd)raw >$(basename $i mhd)raw
done

