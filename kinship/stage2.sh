#!/bin/bash
#$ -S /bin/bash
#$ -o /dev/null
#$ -e /dev/null
#$ -cwd
#$ -V
#$ -l tmem=6G,h_vmem=6G
#$ -l h_rt=24:0:0
#$ -t 1-25
set -u
set -x

scriptname=`basename $0`
mkdir -p ${scriptname}.qsub.out ${scriptname}.qsub.err
exec >${scriptname}.qsub.out/${scriptname}_${SGE_TASK_ID}_${JOB_ID}.out 2>${scriptname}.qsub.err/${scriptname}_${SGE_TASK_ID}_${JOB_ID}.err
args=( header 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y )
chr=${args[$SGE_TASK_ID]}

infile=chr${chr}.recode.vcf
outfile=chr${chr}.vcf.gz

filters=''
bcftools view --samples-file samples.txt $filters --output-file $outfile --output-type z $infile
vcftools --gzvcf $outfile --plink --out chr${chr}_plink 


