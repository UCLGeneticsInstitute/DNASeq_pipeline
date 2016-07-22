#$ -S /bin/sh
#$ -l h_vmem=15G
#$ -l tmem=15G
#$ -l h_rt=50:00:0
#$ -V
#$ -R y
#$ -pe smp 1
#$ -cwd
#$ -o /cluster/project8/vyp/cian/scripts/cluster/output/
#$ -e /cluster/project8/vyp/cian/scripts/cluster/error/
/share/apps/R-3.1.0/bin/R CMD BATCH --no-save plot_singleVariant_results.R
