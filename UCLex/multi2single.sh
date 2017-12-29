
cd ${UCLEX}/bgt
zcat ${UCLEX_RELEASE} | cut -f-8 | python for_vep.py > for_vep.vcf

