#!/bin/bash
for file in $(find /SAN/vyplab/UCLex/mainset_July2016/cian/SKAT -name '*_snps' )
do
    echo $file
    awk 'NR == 1 || $10 < min {line = $0; min = $10}END{print line}' $file
done
