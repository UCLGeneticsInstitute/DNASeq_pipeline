#$ -S /bin/bash
#$ -N SV-AgreeTable-STEP3
#$ -l tmem=20G
#$ -l h_vmem=20G
#$ -pe smp 1
#$ -j y
#$ -R y
#$ -cwd
#$ -o /SAN/vyplab/NCMD/b37/SV-MERGED
#$ -e /SAN/vyplab/NCMD/b37/SV-MERGED
#$ -t 1-24

args=( header 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y )
chr=${args[$SGE_TASK_ID]}

R --vanilla < /SAN/vyplab/NCMD/b37/SV-MERGED/SV-AgreementTable-STEP3.r --args $chr 1
R --vanilla < /SAN/vyplab/NCMD/b37/SV-MERGED/SV-AgreementTable-STEP3.r --args $chr 2
R --vanilla < /SAN/vyplab/NCMD/b37/SV-MERGED/SV-AgreementTable-STEP3.r --args $chr 3
R --vanilla < /SAN/vyplab/NCMD/b37/SV-MERGED/SV-AgreementTable-STEP3.r --args $chr 4
R --vanilla < /SAN/vyplab/NCMD/b37/SV-MERGED/SV-AgreementTable-STEP3.r --args $chr 5
R --vanilla < /SAN/vyplab/NCMD/b37/SV-MERGED/SV-AgreementTable-STEP3.r --args $chr 6
