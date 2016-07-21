#!/bin/bash
#$ -S /bin/bash
#$ -o /dev/null
#$ -e /dev/null
#$ -cwd
#$ -V
#$ -l tmem=2.5G,h_vmem=2.5G
#$ -l h_rt=24:0:0
#$ -t 1-25
set -u
set -x

scriptname=`basename $0`
mkdir -p ${scriptname}.qsub.out ${scriptname}.qsub.err
exec >${scriptname}.qsub.out/${scriptname}_${SGE_TASK_ID}_${JOB_ID}.out 2>${scriptname}.qsub.err/${scriptname}_${SGE_TASK_ID}_${JOB_ID}.err
args=( header 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y )
chr=${args[$SGE_TASK_ID]}

vcf=/goon2/scratch2/vyp-scratch2/vincent/GATK/mainset_February2015b/mainset_February2015b_chr${chr}_indels_filtered.vcf.gz

vcftools --gzvcf $vcf \
--chr ${chr} \
--minGQ 10 \
--minDP 5 \
--max-missing 0.8 \
--maf 0.01 \
--max-maf 0.99 \
--recode \
--out chr${chr} > chr${chr}.log

echo $chr
echo Finished
date
