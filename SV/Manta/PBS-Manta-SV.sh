#$ -N Nik-Manta-SV
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
python=/home/ucbtlcu/bin/0-Python/python
picard=/home/ucbtlcu/bin/picard.jar
configManta=/SAN/mottlab/heterosis/bin/6_Manta/bin/configManta.py
refdir=/SAN/vyplab/UKIRDC/reference
bamdir=/SAN/vyplab/UKIRDC/b37
outdir=/SAN/mottlab/Mitochondria/5_Nik/Manta


cd $outdir

for f in $bamdir/*_sorted_unique.bam
do
 x=$f
 y=${x%_sorted_unique.bam}
 z=${y##*/}
 mkdir $z
 cd $z
 $python $configManta --bam $f --referenceFasta $refdir/human_g1k_v37.fasta --runDir $outdir/$z
 $python runWorkflow.py -m local -j 1
 cd ..
 echo $z
done

date

