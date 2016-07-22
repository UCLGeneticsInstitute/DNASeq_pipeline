#$ -N Nik-SoftSearch-SV
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
SoftSearch=/SAN/mottlab/heterosis/bin/4_SoftSearch/src/SoftSearch.pl
SoftSearch_multi=/SAN/mottlab/heterosis/bin/4_SoftSearch/src/SoftSearch.multi.pl
SoftSearch_Filter=/SAN/mottlab/heterosis/bin/4_SoftSearch/src/SoftSearch_Filter.pl
Annotate_SoftSearch=/SAN/mottlab/heterosis/bin/4_SoftSearch/src/Annotate_SoftSearch.pl
perl=/usr/bin/perl
refdir=/SAN/vyplab/UKIRDC/reference
bamdir=/SAN/vyplab/UKIRDC/b37
outdir=/SAN/mottlab/Mitochondria/5_Nik/SoftSearch

cd $outdir

for f in $bamdir/*_sorted_unique.bam
do
 x=$f
 y=${x%_sorted_unique.bam};
 z=${y##*/};
 $perl $SoftSearch -b $f -f $refdir/human_g1k_v37.fasta -o SoftSearch.$z.vcf
 echo $z
done

date

