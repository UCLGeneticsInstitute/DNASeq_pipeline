#$ -N Nik-Wham-SV
#$ -l tmem=20G
#$ -l h_vmem=20G
#$ -S /bin/bash
#$ -pe smp 1
#$ -j y
#$ -R y
#$ -cwd
#$ -e /home/ucbtlcu/qsublog/
#$ -o /home/ucbtlcu/qsublog/

# Kill script if any commands fail
set -e

date

bedtools=/home/ucbtlcu/bin/0-bedtools2/bin/bedtools
bwa=/home/ucbtlcu/bin/0-bwa-0.7.15/bwa
samtools=/share/apps/genomics/samtools-1.2/bin/samtools
bcftools=/home/ucbtlcu/bin/bcftools
java=/share/apps/java/bin/java
picard=/home/ucbtlcu/bin/picard.jar
wham=/SAN/mottlab/heterosis/bin/1-Wham/bin/wham
whamg=/SAN/mottlab/heterosis/bin/1-Wham/bin/whamg
refdir=/SAN/vyplab/UKIRDC/reference
bamdir=/SAN/vyplab/UKIRDC/b37
outdir=/SAN/mottlab/Mitochondria/5_Nik/Wham

cd $outdir

for f in $bamdir/*_sorted_unique.bam
do
 x=$f
 y=${x%_sorted_unique.bam};
 z=${y##*/};
 $wham -x 1 -f $refdir/human_g1k_v37.fasta -t $f > WHAM.$z.vcf 2> WHAM.$z.err
 echo $z
done

ls $bamdir/*_sorted_unique.bam > $outdir/bamlist
$whamg -g WHAMG.graph.txt -x 1 -a $refdir/human_g1k_v37.fasta -c 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y -f bamlist > WHAMG.vcf 2> WHAMG.err

date

