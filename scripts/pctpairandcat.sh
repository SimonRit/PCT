#!/bin/sh -x

for FILENAME in results.???? 
do
    cd $FILENAME
    $(dirname $0)/pctpairprotons.sh
    cd ..
done

for FILENAME in $(ls */pairs0000.mhd)
do
    if ! test -e merged
    then
        mkdir merged
        cp $(dirname $FILENAME)/pairs????.mhd merged
    else
        for i in $(dirname $FILENAME)/pairs????.mhd
        do
            CURRENTSIZE=$(grep DimSize merged/$(basename $i) | sed 's/.*DimSize = \(.*\) \(.*\)/\2/g' )
            NEWSIZE=$(grep DimSize $i | sed 's/.*DimSize = \(.*\) \(.*\)/\2/g' )
            SUMSIZE=$(($CURRENTSIZE+$NEWSIZE))
            sed -i "s/.*DimSize = \(.*\) \(.*\)/DimSize = \1 $(($CURRENTSIZE+$NEWSIZE))/g" merged/$(basename $i)
        done
    fi
done

cd merged
for i in pairs????.mhd
do
    cat ../results.????/$(basename $i mhd)raw >$(basename $i mhd)raw
done

