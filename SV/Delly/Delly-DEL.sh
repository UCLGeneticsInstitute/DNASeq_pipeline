#$ -N Nik-Delly-DEL
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
delly=/SAN/mottlab/heterosis/bin/2_Delly/pre-compiled/delly_v0.7.3_linux_x86_64bit
delly_parallel=/SAN/mottlab/heterosis/bin/2_Delly/pre-compiled/delly_v0.7.3_parallel_linux_x86_64bit
refdir=/SAN/vyplab/UKIRDC/reference
bamdir=/SAN/vyplab/UKIRDC/b37
outdir=/SAN/mottlab/Mitochondria/5_Nik/Delly/DEL

cd $outdir

for f in $bamdir/*_sorted_unique.bam
do
 x=$f
 y=${x%_sorted_unique.bam}
 z=${y##*/}
 $delly call -t DEL -g $refdir/human_g1k_v37.fasta $f -o $z.bcf
 echo $z
done

$delly merge -t DEL -m 0 -n 1000000 -o DEL1.bcf -b 1000 -r 0.8 chr5_FR07919256.bcf chr5_FR07922004.bcf

for f in $bamdir/*.bam
do
 x=$f
 y=${x%_sorted_unique.bam}
 z=${y##*/}
 $delly call -t DEL -g $refdir/human_g1k_v37.fasta $f -v DEL1.bcf -o $z.geno.bcf
 echo $z
done

$bcftools merge -O b -o DEL2.bcf chr5_FR07919256.geno.bcf chr5_FR07922004.geno.bcf
$bcftools index DEL2.bcf DEL2
$delly filter -t DEL -f germline -o DEL3.bcf -g $refdir/human_g1k_v37.fasta DEL2.bcf

date

