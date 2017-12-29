
cd bgt
bgt import -S ${UCLEX_RELEASE}.bgt ${UCLEX_RELEASE}.vcf.gz 

cp ${UCLEX_RELEASE}.bgt.spl ${UCLEX_RELEASE}.bgt.spl.bak
awk -v FS="\t" -v OFS="\t" '{print $0,"name:Z:"$1}' ${UCLEX_RELEASE}.bgt.spl.bak > ${UCLEX_RELEASE}.bgt.spl


