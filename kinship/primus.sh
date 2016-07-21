#$ -S /bin/bash
#$ -o /dev/null
#$ -e /dev/null
#$ -cwd
#$ -V
#$ -l tmem=8G,h_vmem=8G
#$ -l h_rt=24:0:0

set -u
set -x
scriptname=primus
mkdir -p ${scriptname}.qsub.out ${scriptname}.qsub.err
exec >${scriptname}.qsub.out/${scriptname}_${SGE_TASK_ID}_${JOB_ID}.out 2>${scriptname}.qsub.err/${scriptname}_${SGE_TASK_ID}_${JOB_ID}.err


#plink --noweb --file chr1_plink --merge-list chromosomes_2-22.txt --make-bed --out all 
#plink --noweb --bfile all --geno 0.1 --mind 0.1 --make-bed --out all_filtered

king -b all.bed --kinship --prefix all


tail -n+2 all.kin0 | awk '$1!=$3' | sort -k8nr  > all.kin0_sorted2 
cat <(head -n1 all.kin0) all.kin0_sorted2 > all.kin0_sorted

#grep -v Levine king_all_filtered.kin0 > king_all_filtered.kin0_no_Levine

#PRIMUS=/cluster/project8/vyp/pontikos/Software/PRIMUS_v1.8.0/bin/run_PRIMUS.pl
#perl $PRIMUS -i FILE=king_all_filtered.kin0_sorted -threshold 0.1

PRIMUS=/cluster/project8/vyp/AdamLevine/software/PRIMUS/PRIMUSv0.5/bin/PRIMUS-v0.5.pl
perl $PRIMUS -input all.kin0_sorted -threshold 0.1

independent=all.kin0_sorted_results/all.kin0_sorted_maximum_independent_set
grep -v OneKG $independent > independent_noOneKG.txt
awk 'NR>1 {print $1}' independent_noOneKG.txt > UCL-exome_unrelated.txt

