mkdir in out 2>/dev/null
DIMX=$(grep DimSize proj000.mhd | sed "s/DimSize = \(.*\) \(.*\) .*/\1/g")
DIMY=$(grep DimSize proj000.mhd | sed "s/DimSize = \(.*\) \(.*\) .*/\2/g")
ORIGINZ=$(grep Offset proj000.mhd | sed "s/Offset = \(.*\) \(.*\) \(.*\)/\3/g")
SPACINGZ=$(grep ElementSpacing proj000.mhd | sed "s/ElementSpacing = \(.*\) \(.*\) \(.*\)/\3/g")
SLICEIN=$(echo "scale=0 ; (-100 - ${ORIGINZ})/${SPACINGZ}" | bc)
SLICEOUT=$(echo "scale=0 ; (100 - ${ORIGINZ})/${SPACINGZ}" | bc)
for i in proj???.mhd
do
    echo HeaderSize = $((4*${DIMX}*${DIMY}*${SLICEIN}))  > in/${i}
    echo HeaderSize = $((4*${DIMX}*${DIMY}*${SLICEOUT})) > out/${i}

    for f in in out
    do
    cat $i | sed "s/proj/\.\.\/proj/g" \
           | sed "s/NDims = 3/NDims = 2/g" \
           | sed "s/CenterOfRotation.*/CenterOfRotation = 0 0/g" \
           | sed "s/TransformMatrix.*/TransformMatrix = 1 0 0 1/g" \
           | sed "s/DimSize = \(.*\) \(.*\) .*/DimSize = \1 \2/g" \
           | sed "s/ElementSpacing = \(.*\) \(.*\) .*/ElementSpacing = \1 \2/g" \
           | sed "s/Offset = \(.*\) \(.*\) .*/Offset = \1 \2/g" \
           >>${f}/${i}
    done
done

