#$ -N Nik-Lumpy-SV
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
lumpy=/SAN/mottlab/heterosis/bin/3_Lumpy/bin/lumpy
lumpyexpress=/SAN/mottlab/heterosis/bin/3_Lumpy/bin/lumpyexpress
extractSplitReads_BwaMem=/SAN/mottlab/heterosis/bin/3_Lumpy/scripts/extractSplitReads_BwaMem
svtyper=/home/ucbtlcu/bin/0-svtyper-master/svtyper
refdir=/SAN/vyplab/UKIRDC/reference
bamdir=/SAN/vyplab/UKIRDC/b37
outdir=/SAN/mottlab/Mitochondria/5_Nik/Lumpy

cd $outdir

for f in $bamdir/*_sorted_unique.bam
do
 x=$f
 y=${x%_sorted_unique.bam}
 z=${y##*/}
 $samtools view -b -F 1294 $f > tmp.discordants.bam
 $samtools view -h $f | $extractSplitReads_BwaMem -i stdin | $samtools view -Sb - > tmp.splitters.bam
 $samtools sort tmp.discordants.bam tmp.discordants.sorted
 $samtools sort tmp.splitters.bam tmp.splitters.sorted
 $samtools index tmp.splitters.sorted.bam
 $lumpyexpress -B $f -S tmp.splitters.sorted.bam -D tmp.discordants.sorted.bam -o $z.vcf
 #$svtyper -B $f -S tmp.splitters.sorted.bam -i $z.vcf > $z.gt.vcf
 rm tmp.*
 echo $z
done

date

