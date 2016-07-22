#$ -N Nik-SpeedSeq-SV
#$ -l tmem=10G
#$ -l h_vmem=10G
#$ -S /bin/bash
#$ -pe smp 2
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
speedseq=/SAN/mottlab/heterosis/bin/7_SpeedSeq/bin/speedseq
refdir=/SAN/vyplab/UKIRDC/reference
bamdir=/SAN/vyplab/UKIRDC/b37
outdir=/SAN/mottlab/Mitochondria/5_Nik/SpeedSeq


cd $outdir

for f in $bamdir/*_sorted_unique.bam
do
 x=$f
 y=${x%_sorted_unique.bam};
 z=${y##*/};
 $speedseq var -t 2 -o VAR.$z $refdir/human_g1k_v37.fasta $f
 echo $z
done

for f in $bamdir/*_sorted_unique.bam
do
 x=$f
 y=${x%_sorted_unique.bam};
 z=${y##*/};
 $speedseq sv -t 2 -P -o SV.$z -B $f -R $refdir/human_g1k_v37.fasta
 echo $z
done

date

