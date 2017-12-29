# split single lines
zcat ${output}_chr${chr}_for_annovar.vcf.gz | /share/apps/python/bin/python ${baseFolder}/annotation/multiallele_to_single_vcf.py --headers CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO > ${output}_VEP/chr${chr}_for_VEP.vcf
/share/apps/genomics/htslib-1.1/bin/bgzip -f -c ${output}_VEP/chr${chr}_for_VEP.vcf > ${output}_VEP/chr${chr}_for_VEP.vcf.gz
/share/apps/genomics/htslib-1.1/bin/tabix -f -p vcf ${output}_VEP/chr${chr}_for_VEP.vcf.gz
rm ${output}_VEP/chr${chr}_for_VEP.vcf
