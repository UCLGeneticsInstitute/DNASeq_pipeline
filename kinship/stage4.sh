#!/bin/bash

#$ -S /bin/bash
#$ -o /home/zchads1/cluster/UCL-exomes/log/
#$ -e /home/zchads1/cluster/UCL-exomes/log/
#$ -l h_rt=2:0:0
#$ -l tmem=6G,h_vmem=6G
#$ -t 1-23

i=$(($SGE_TASK_ID - 1))

date

cd /home/zchads1/cluster/UCL-exomes/freq

UNRELATED=../common_good_plink/UCL-exome_unrelated.txt

chromosomes=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X)
chr=${chromosomes[$i]}

vcftools=~/cluster/software/vcftools_0.1.12b/bin/vcftools

vcf=../by_chr/chr_${chr}.vcf.gz

$vcftools --gzvcf $vcf \
--keep $UNRELATED \
--freq \
--out chr${chr} > chr${chr}.log


echo $chr
echo Finished
date

