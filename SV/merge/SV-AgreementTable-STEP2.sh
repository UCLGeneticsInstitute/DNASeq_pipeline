#$ -S /bin/bash
#$ -N SV-AgreeTable-STEP2
#$ -l tmem=30G
#$ -l h_vmem=30G
#$ -pe smp 1
#$ -j y
#$ -R y
#$ -cwd
#$ -o /SAN/vyplab/NCMD/b37/SV-MERGED
#$ -e /SAN/vyplab/NCMD/b37/SV-MERGED
#$ -t 1-6

R --vanilla < /SAN/vyplab/NCMD/b37/SV-MERGED/SV-AgreementTable-STEP2.r --args $SGE_TASK_ID
