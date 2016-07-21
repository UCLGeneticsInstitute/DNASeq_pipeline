#!/bin/bash
#$ -S /bin/bash
#$ -o /dev/null
#$ -e /dev/null
#$ -cwd
#$ -V
#$ -l scr=1G
#$ -l tmem=2G,h_vmem=2G
#$ -l h_rt=2:0:0
#$ -t 1-25
set -u
set -x
scriptname=king
mkdir -p ${scriptname}.qsub.out ${scriptname}.qsub.err
exec >${scriptname}.qsub.out/${scriptname}_${SGE_TASK_ID}_${JOB_ID}.out 2>${scriptname}.qsub.err/${scriptname}_${SGE_TASK_ID}_${JOB_ID}.err
args=( header 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y )
f=${args[$SGE_TASK_ID]}


king -b chr${f}.bed --kinship --related --degree 2 --prefix chr${f}

