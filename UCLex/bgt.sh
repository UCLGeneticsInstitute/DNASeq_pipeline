
bcftools concat ${UCLEX_RELEASE}_chr1.vcf.gz ${UCLEX_RELEASE}_chr2.vcf.gz ${UCLEX_RELEASE}_chr3.vcf.gz ${UCLEX_RELEASE}_chr4.vcf.gz ${UCLEX_RELEASE}_chr5.vcf.gz ${UCLEX_RELEASE}_chr6.vcf.gz ${UCLEX_RELEASE}_chr7.vcf.gz ${UCLEX_RELEASE}_chr8.vcf.gz ${UCLEX_RELEASE}_chr9.vcf.gz ${UCLEX_RELEASE}_chr10.vcf.gz ${UCLEX_RELEASE}_chr11.vcf.gz ${UCLEX_RELEASE}_chr12.vcf.gz ${UCLEX_RELEASE}_chr13.vcf.gz ${UCLEX_RELEASE}_chr14.vcf.gz ${UCLEX_RELEASE}_chr15.vcf.gz ${UCLEX_RELEASE}_chr16.vcf.gz ${UCLEX_RELEASE}_chr17.vcf.gz ${UCLEX_RELEASE}_chr18.vcf.gz ${UCLEX_RELEASE}_chr19.vcf.gz ${UCLEX_RELEASE}_chr20.vcf.gz ${UCLEX_RELEASE}_chr21.vcf.gz ${UCLEX_RELEASE}_chr22.vcf.gz ${UCLEX_RELEASE}_chrX.vcf.gz -O z -o bgt/${UCLEX_RELEASE}.vcf.gz 

cd bgt
bgt import -S ${UCLEX_RELEASE}.bgt ${UCLEX_RELEASE}.vcf.gz 

cp ${UCLEX_RELEASE}.bgt.spl ${UCLEX_RELEASE}.bgt.spl.bak
awk -v FS="\t" -v OFS="\t" '{print $0,"name:Z:"$1}' ${UCLEX_RELEASE}.bgt.spl.bak > ${UCLEX_RELEASE}.bgt.spl


