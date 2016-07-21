#!/bin/bash

#$ -S /bin/bash
#$ -o /home/zchads1/cluster/UCL-exomes/log/
#$ -e /home/zchads1/cluster/UCL-exomes/log/
#$ -l h_rt=1:0:0
#$ -l tmem=2G,h_vmem=2G
#$ -t 1-23

i=$(($SGE_TASK_ID - 1))

date

cd /home/zchads1/cluster/UCL-exomes/freq

chromosomes=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X)
chr=${chromosomes[$i]}

bgzip=/cluster/project8/vyp/vincent/Software/tabix-0.2.5/bgzip
tabix=/cluster/project8/vyp/vincent/Software/tabix-0.2.5/tabix

python ../prepare_af_from_vcftools.py chr${chr}.frq  | $bgzip > UCLfreq_${chr}.vcf.gz
$tabix -p vcf UCLfreq_${chr}.vcf.gz

echo $chr
echo Finished
date

