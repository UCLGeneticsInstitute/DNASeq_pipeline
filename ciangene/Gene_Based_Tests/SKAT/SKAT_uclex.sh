#$ -S /bin/sh
#$ -l h_vmem=10G
#$ -l tmem=10G
#$ -l h_rt=2:00:0
#$ -V
#$ -R y
#$ -pe smp 1
#$ -cwd
#$ -o /cluster/project8/vyp/cian/scripts/cluster/output/
#$ -e /cluster/project8/vyp/cian/scripts/cluster/error/
/share/apps/R/bin/R CMD BATCH --no-save SKAT_uclex.R
