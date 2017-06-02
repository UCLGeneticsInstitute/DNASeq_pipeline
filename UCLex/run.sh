set -x
#set -e

chr=$1
release=$2

memoSmall=10

# prints to stderr in red
function error() { >&2 echo -e "\033[31m$*\033[0m"; }
function stop() { error "$*"; exit 1; }
try() { "$@" || stop "cannot $*"; }

mainFolder=/SAN/vyplab/UCLex

referenceFolder=/cluster/scratch3/vyp-scratch2
#fasta=${referenceFolder}/reference_datasets/human_reference_sequence/human_g1k_v37.fasta
fasta=/SAN/vyplab/UKIRDC/reference/human_g1k_v37.fasta
bundle=${referenceFolder}/reference_datasets/GATK_bundle
#dbsnp=${bundle}/dbsnp_137.b37.vcf
dbsnp=/SAN/vyplab/UKIRDC/reference/dbsnp_138.b37.vcf.gz
#Rscript=/share/apps/R-3.3.0/bin/Rscript
#Rbin=/share/apps/R-3.3.0/bin/R 
Rbin=/cluster/project8/vyp/vincent/Software/R-3.3.0/bin/R
Rscript=/cluster/project8/vyp/vincent/Software/R-3.3.0/bin/Rscript

#java=/share/apps/jdk1.7.0_45/bin/java
java=/share/apps/jdk/jre/bin/java
tmpDir=/scratch0/vyp
target=/cluster/project8/vyp/exome_sequencing_multisamples/target_region/data/merged_exome_target_cleaned.bed

GATK=/cluster/project8/vyp/vincent/Software/GenomeAnalysisTK-3.5-0/GenomeAnalysisTK.jar

baseFolder="$( cd "$( dirname "${BASH_SOURCE[0]}" )/.." && pwd )"
echo $baseFolder

gVCFlist=gvcf_list.mainset_${release}

maxGaussians=6
maxGaussiansIndels=5
numBad=1000
numBadIndels=1000
GQ=20


output=${mainFolder}/mainset_${release}/mainset_${release}

if [ ! -s mainset_${release}_gVCF_${chr}.list ]
then
    echo "Make gVCF file list"
    while read path id format
        do
            gVCF=${path}/chr${chr}/${id}
            if [ ! -s $gVCF ]; then stop "Cannot find $gVCF"; fi
            if [ ! -s $gVCF.tbi ]; then stop "Cannot find $gVCF.tbi"; fi
            echo "$gVCF" 
    done < <(tail -n +2 $gVCFlist) > mainset_${release}_gVCF_${chr}.list
fi

echo "Running the genotype module"
script=${scripts_folder}/subscript_chr${chr}.sh
f=${output}_chr${chr}.vcf.gz.tbi

if [ ! -s  ${output}_chr${chr}.vcf.gz ]
then
    echo GenotypeGVCFs
    $java -Djava.io.tmpdir=/scratch0/ -Xmx${memoSmall}g -jar $GATK -R $fasta -L $chr -L $target -T GenotypeGVCFs \
    --interval_set_rule INTERSECTION --interval_padding 100  \
    --annotation InbreedingCoeff --annotation QualByDepth --annotation HaplotypeScore \
    --annotation MappingQualityRankSumTest --annotation ReadPosRankSumTest --annotation FisherStrand \
    --dbsnp ${dbsnp} -V mainset_${release}_gVCF_${chr}.list --out ${output}_chr${chr}.vcf.gz 
fi

if [ ! -s ${output}_chr${chr}_SNPs.vcf.gz ]
then
    echo extract the SNPs: SelectVariants
    $java  -Djava.io.tmpdir=/scratch0/ -Xmx${memoSmall}g -jar ${GATK} -R $fasta -L $chr -T SelectVariants \
-selectType SNP -V ${output}_chr${chr}.vcf.gz  --out ${output}_chr${chr}_SNPs.vcf.gz
fi

if [ ! -s ${output}_chr${chr}.vcf.gz --out ${output}_chr${chr}_indels.vcf.gz ]
then
    echo extract the indels
    $java  -Djava.io.tmpdir=/scratch0/ -Xmx${memoSmall}g -jar ${GATK}  -R $fasta  -L $chr  -T SelectVariants \
    -selectType INDEL \
    -selectType MIXED \
    -V ${output}_chr${chr}.vcf.gz --out ${output}_chr${chr}_indels.vcf.gz
fi

echo VariantRecalibrator
echo SNP recalibration
for maxGauLoc in $(seq 3 6 | sort -r)
do
    if [[ -s ${output}_chr${chr}_SNPs_filtered.vcf.gz ]]
    then
        break
    fi
    $java -Djava.io.tmpdir=/scratch0/ -Xmx${memoSmall}g -jar ${GATK} -R $fasta  -L $chr -T VariantRecalibrator \
       --input ${output}_chr${chr}_SNPs.vcf.gz --maxGaussians ${maxGauLoc} --mode SNP \
       -resource:hapmap,VCF,known=false,training=true,truth=true,prior=15.0 ${bundle}/hapmap_3.3.b37.vcf  \
       -resource:omni,VCF,known=false,training=true,truth=false,prior=12.0 ${bundle}/1000G_omni2.5.b37.vcf \
       -resource:dbsnp,VCF,known=true,training=false,truth=false,prior=8.0 ${bundle}/dbsnp_137.b37.vcf \
       -an QD -an FS -an ReadPosRankSum -an InbreedingCoeff \
       -tranche 100.0 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 \
       --minNumBadVariants ${numBad} \
       -recalFile ${output}_chr${chr}_SNPs_combrec.recal \
       -tranchesFile ${output}_chr${chr}_SNPs_combtranch \
       -rscriptFile  ${output}_chr${chr}_recal_plots_snps.R

    /share/apps/R/bin/Rscript ${output}_chr${chr}_recal_plots_snps.R
    # apply_recal: ApplyRecalibration
    $java -Xmx${memoSmall}g -jar ${GATK} -R $fasta -L $chr  -T ApplyRecalibration \
       --ts_filter_level 99.5 \
       --recal_file ${output}_chr${chr}_SNPs_combrec.recal \
       --tranches_file ${output}_chr${chr}_SNPs_combtranch --mode SNP \
       --input ${output}_chr${chr}_SNPs.vcf.gz --out ${output}_chr${chr}_SNPs_filtered.vcf.gz
done


if [ ! -s ${output}_chr${chr}_indels_hard_filtered.vcf.gz.tbi ]
then
    echo Indel hard filtering
    $java -Djava.io.tmpdir=/scratch0/ -Xmx${memoSmall}g -jar ${GATK} -R $fasta  -L $chr  -T VariantFiltration \
   --filterExpression "QD < 2.0 || FS > 50.0 || ReadPosRankSum < -20.0" \
   --filterName "FAIL" \
   -V ${output}_chr${chr}_indels.vcf.gz --out ${output}_chr${chr}_indels_hard_filtered.vcf.gz
fi

if [ ! -s ${output}_chr${chr}_indels_filtered.vcf.gz ]
then
    echo Indel recalibration
    # calculate recal: VariantRecalibrator
    for maxGaussiansIndels in $(seq 3 6 | sort -r)
    do
        if [[ -s ${output}_chr${chr}_indels_filtered.vcf.gz ]]
        then
            break
        fi
    $java -Djava.io.tmpdir=${tmpDir} -Xmx${memoSmall}g -jar ${GATK} -R $fasta  -L $chr -T VariantRecalibrator \
       -resource:mills,known=true,training=true,truth=true,prior=12.0 ${bundle}/Mills_and_1000G_gold_standard.indels.b37.vcf \
       -an QD -an FS -an ReadPosRankSum -an InbreedingCoeff \
       -tranche 100.0 -tranche 99.5  -tranche 99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \
       --minNumBadVariants ${numBadIndels} \
       --maxGaussians ${maxGaussiansIndels} \
       --input ${output}_chr${chr}_indels.vcf.gz  --mode INDEL \
       -recalFile ${output}_chr${chr}_indels_combrec.recal \
       -tranchesFile ${output}_chr${chr}_indels_combtranch \
       -rscriptFile ${output}_chr${chr}_recal_plots_indels.R 
    /share/apps/R/bin/Rscript ${output}_chr${chr}_recal_plots_indels.R
    #apply_recal: ApplyRecalibration
    $java -Djava.io.tmpdir=${tmpDir} -Xmx${memoSmall}g -jar ${GATK} -R $fasta  -L $chr -T ApplyRecalibration  \
       --ts_filter_level 95.0 \
       --recal_file ${output}_chr${chr}_indels_combrec.recal \
       --tranches_file ${output}_chr${chr}_indels_combtranch --mode INDEL \
       --input ${output}_chr${chr}_indels.vcf.gz --out ${output}_chr${chr}_indels_filtered.vcf.gz
    done
fi


if [ ! -s ${output}_chr${chr}_filtered.vcf.gz.tbi ]
then
    echo Recal merge SNPs and indels
    $java -Djava.io.tmpdir=${tmpDir} -Xmx${memoSmall}g -jar ${GATK} -R $fasta -L $chr -T CombineVariants \
   --assumeIdenticalSamples --variant:SNPs ${output}_chr${chr}_SNPs_filtered.vcf.gz \
   --variant:indels ${output}_chr${chr}_indels_filtered.vcf.gz \
   -genotypeMergeOptions PRIORITIZE  \
   -priority SNPs,indels \
   --out ${output}_chr${chr}_filtered.vcf.gz
fi


if [ ! -s ${output}_chr${chr}_for_annovar.vcf.gz.tbi ]
then
    echo annovar chr${chr}
    echo creates ${output}_chr${chr}_recal_filtered2.vcf 
    zcat ${output}_chr${chr}_filtered.vcf.gz | perl ${baseFolder}/UCLex/custom_filtering.pl ${output}_chr${chr}_filtered2.vcf ${GQ}
    cut -f1-8 ${output}_chr${chr}_filtered2.vcf > ${output}_chr${chr}_for_annovar.vcf
    /cluster/project8/vyp/vincent/Software_heavy/annovar_Feb2013/convert2annovar.pl --allallele -format vcf4 --includeinfo ${output}_chr${chr}_for_annovar.vcf > ${output}_chr${chr}_db
    /cluster/project8/vyp/vincent/Software_heavy/annovar_Feb2013/summarize_annovar_VP.pl -ver1000g 1000g2012apr -verdbsnp 137 -veresp 6500si -alltranscript -buildver hg19 --genetype ensgene --remove ${output}_chr${chr}_db /cluster/project8/vyp/vincent/Software_heavy/annovar_Feb2013/humandb_hg19/ 
    python ${baseFolder}/UCLex/annovar_vcf_combine_VP.py ${output}_chr${chr}_filtered2.vcf ${output}_chr${chr}_db.exome_summary.csv ${output}_chr${chr}_exome_table.csv
    perl ${baseFolder}/UCLex/make_matrix_calls.pl ${output}_chr${chr}_exome_table.csv ${output} $chr
    echo compress ${output}_chr${chr}_filtered2.vcf 
    /share/apps/genomics/htslib-1.1/bin/bgzip -f -c ${output}_chr${chr}_filtered2.vcf > ${output}_chr${chr}_filtered2.vcf.gz
    /share/apps/genomics/htslib-1.1/bin/tabix -f -p vcf ${output}_chr${chr}_filtered2.vcf.gz
    rm ${output}_chr${chr}_filtered2.vcf
    echo compress ${output}_chr${chr}_for_annovar.vcf
    /share/apps/genomics/htslib-1.1/bin/bgzip -f -c ${output}_chr${chr}_for_annovar.vcf > ${output}_chr${chr}_for_annovar.vcf.gz
    /share/apps/genomics/htslib-1.1/bin/tabix -f -p vcf ${output}_chr${chr}_for_annovar.vcf.gz
    rm ${output}_chr${chr}_for_annovar.vcf
fi


if [ ! -s ${output}_VEP/chr${chr}_for_VEP.vcf.gz.tbi ]
then
    echo VEP_input
    # split single lines
    zcat ${output}_chr${chr}_for_annovar.vcf.gz | /share/apps/python/bin/python ${baseFolder}/annotation/multiallele_to_single_vcf.py --headers CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO > ${output}_VEP/chr${chr}_for_VEP.vcf
    /share/apps/genomics/htslib-1.1/bin/bgzip -f -c ${output}_VEP/chr${chr}_for_VEP.vcf > ${output}_VEP/chr${chr}_for_VEP.vcf.gz
    /share/apps/genomics/htslib-1.1/bin/tabix -f -p vcf ${output}_VEP/chr${chr}_for_VEP.vcf.gz
    rm ${output}_VEP/chr${chr}_for_VEP.vcf
fi

mkdir -p ${output}_VEP


if [ ! -s ${output}_VEP/CADD_chr${chr}.vcf.gz.tbi ]
then
    echo  CADD
    # split single lines
    zcat ${output}_chr${chr}_for_annovar.vcf.gz | /share/apps/python/bin/python ${baseFolder}/annotation/multiallele_to_single_vcf.py --headers CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO > ${output}_VEP/chr${chr}_for_VEP.vcf
    /share/apps/genomics/htslib-1.1/bin/bgzip -f -c ${output}_VEP/chr${chr}_for_VEP.vcf > ${output}_VEP/chr${chr}_for_VEP.vcf.gz
    /share/apps/genomics/htslib-1.1/bin/tabix -f -p vcf ${output}_VEP/chr${chr}_for_VEP.vcf.gz
    rm ${output}_VEP/chr${chr}_for_VEP.vcf
    echo  RUN CADD
    /share/apps/genomics/CADD_v1.3/bin/score.sh ${output}_VEP/chr${chr}_for_VEP.vcf.gz ${output}_VEP/CADD_chr${chr}.vcf.gz
    /share/apps/genomics/htslib-1.1/bin/tabix -p vcf ${output}_VEP/CADD_chr${chr}.vcf.gz
fi

mkdir -p ${output}_snpStats
echo convertToR chr${chr}
/share/apps/R/bin/Rscript ${baseFolder}/UCLex/convert_to_R.R --chromosome ${chr} --root ${output} > ${output}_snpStats/convert_to_R_chr${chr}.out

echo finalCrunch
keyWords=data/controlKeywords.tab
casekeyWords=data/caseKeywords.tab
script=cluster/submission/subscript_chr${chr}.sh
/share/apps/perl/bin/perl ${baseFolder}/UCLex/crunch_controls.pl ${output}_chr${chr}_exome_table.csv $keyWords $casekeyWords ${output}_chr${chr}_exome_crunched.csv data/sampleList_exome.tab none no  ##include all samples


