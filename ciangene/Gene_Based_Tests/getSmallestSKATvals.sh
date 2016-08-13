#!/bin/bash
for file in $(find /SAN/vyplab/UCLex/mainset_July2016/cian/SKAT2 -name '*_snps' )
do
    echo $file
    awk 'NR == 1 || $5 < min {line = $0; min = $5}END{print line}' $file
done
