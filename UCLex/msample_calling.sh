
# prints to stderr in red
function error() { >&2 echo -e "\033[31m$*\033[0m"; }
function stop() { error "$*"; exit 1; }
try() { "$@" || stop "cannot $*"; }

mainFolder=/SAN/vyplab/UCLex
referenceFolder=/cluster/scratch3/vyp-scratch2

fasta=${referenceFolder}/reference_datasets/human_reference_sequence/human_g1k_v37.fasta
bundle=${referenceFolder}/reference_datasets/GATK_bundle

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

submit=no 
force=no
genotype=no
Recal=no
gVCFlist=none

maxGaussians=6
maxGaussiansIndels=5
numBad=1000
numBadIndels=1000
GQ=20

until [ -z "$1" ]
do
    # use a case statement to test vars. we always test $1 and shift at the end of the for block.
    case $1 in
    --force )
        shift
        force=$1;;
    --target )
       shift
       target=$1;;
    --mode )
        shift
        mode=$1;;
    --gVCFlist )
        shift
        gVCFlist=$1;;
    --currentUCLex )
        shift
        currentUCLex=$1;;
    -* )
       echo "Unrecognized option: $1"
       exit 1;;
    esac
    shift
    if [ "$#" = "0" ]; then break; fi
done


############## Now options are all set
scripts_folder=${mainFolder}/mainset_${currentUCLex}/scripts/
stdout_folder=${mainFolder}/mainset_${currentUCLex}/out/
stderr_folder=${mainFolder}/mainset_${currentUCLex}/err/
mkdir -p ${scripts_folder}
mkdir -p ${stdout_folder}
mkdir -p ${stderr_folder}

output=${mainFolder}/mainset_${currentUCLex}/mainset_${currentUCLex}

### Check format of support file.
##should accept tab or space as delimiters
## but does read support tabs and delimeters?
mustBePath=`head -n1 $gVCFlist | cut -f1 -d' ' | cut -f1`
mustBeId=`head -n1 $gVCFlist | cut -f2 -d' ' | cut -f2`

if [[ "$mustBePath" != "path" ]]; then stop "The first column of the file $gVCFlist must have the name path $mustBePath"; fi
if [[ "$mustBeId" != "id" ]]; then stop "The second column of the file $gVCFlist must have the name id $mustBeId"; fi

# set memory
memoSmall=10
memo=15
if [[ "$convertToR" == "yes" || "$annovar" == "yes" ]]; then memo=21.9; fi

mainScript=${scripts_folder}/calling.sh
## individual scripts of the form cluster/submission/subscript_chr${chr}.sh

######## first clean up the individual files
for chr in `seq 1 22` X
do
    f=${scripts_folder}/subscript_chr${chr}.sh
    if [ -e $f ]
    then
        rm $f
    fi
done

#rm cluster/error/*
#rm cluster/out/*

##################################################
function genotype() {
    echo "Running the genotype module"
    for chr in `seq 1 22` X
    do
    script=${scripts_folder}/subscript_chr${chr}.sh
    f=${output}_chr${chr}.vcf.gz.tbi
    if [[ ! -e $f || "$force"=="yes" ]]
    then 
echo "
############### genotype chr${chr}
$java -Djava.io.tmpdir=/scratch0/ -Xmx${memoSmall}g -jar $GATK \\
   -R $fasta \\
   -T GenotypeGVCFs \\
   -L $chr -L $target --interval_set_rule INTERSECTION --interval_padding 100  \\
   --annotation InbreedingCoeff --annotation QualByDepth --annotation HaplotypeScore \\
   --annotation MappingQualityRankSumTest --annotation ReadPosRankSumTest --annotation FisherStrand \\
   --dbsnp ${bundle}/dbsnp_137.b37.vcf \\" >> $script
    while read path id format
        do
            if [[ "$format" == "v1" ]]; then gVCF=${path}/chr${chr}/${id}.gvcf.gz; fi
            if [[ "$format" == "v2" ]]; then gVCF=${path}/${id}-chr${chr}.gvcf.gz; fi
            if [[ "$format" == "v3" ]]; then gVCF=${path}/chr${chr}.gvcf.gz; fi
            if [[ "$format" == "v4" ]]; then gVCF=${path}/chr${chr}/${id}; fi
            echo "Including $gVCF"
            if [ ! -s $gVCF ]; then stop "Cannot find $gVCF"; fi
            if [ ! -s $gVCF.tbi ]; then stop "Cannot find $gVCF.tbi"; fi
            echo "   --variant $gVCF \\" >> $script
        done < <(tail -n +2 $gVCFlist)
echo "   -o ${output}_chr${chr}.vcf.gz" >> $script
    else
        echo $f exists
    fi
    done 
}

##################################################
function recal_extract() {
    for chr in `seq 1 22` X
    do
        #### creates the tmpDir if needed
        tmpDir=/scratch0/GATK_chr${chr}
        f1=${output}_chr${chr}_SNPs.vcf.gz
	f2=${output}_chr${chr}_indels.vcf.gz
	rm -f ${scripts_folder}/subscript_chr${chr}.sh
	if [[ ! -s $f1 ]] || [[ ! -s $f2 ]]
	then
	echo "
############### extract chr${chr}
set +x
mkdir -p $tmpDir
####### first SNPs
#### extract the SNPs: SelectVariants
$java  -Djava.io.tmpdir=${tmpDir} -Xmx${memoSmall}g -jar ${GATK} -R $fasta -L $chr \
   -T SelectVariants \
   -selectType SNP \
   -V ${output}_chr${chr}.vcf.gz  --out ${output}_chr${chr}_SNPs.vcf.gz
####### now indels
#### extract the indels
$java  -Djava.io.tmpdir=${tmpDir} -Xmx${memoSmall}g -jar ${GATK}  -R $fasta  -L $chr \
   -T SelectVariants \
   -selectType INDEL \
   -selectType MIXED \
   -V ${output}_chr${chr}.vcf.gz --out ${output}_chr${chr}_indels.vcf.gz" > ${scripts_folder}/subscript_chr${chr}.sh
	fi
    done
}

###################################################
function recal_snps() {
for chr in `seq 1 22` X
    do
        #### creates the tmpDir if needed
        tmpDir=/scratch0/GATK_chr${chr}
        rm -f ${scripts_folder}/subscript_chr${chr}.sh
	f=${output}_chr${chr}_SNPs_filtered.vcf.gz
	if [[ ! -s $f ]]
        then
	echo "
############### recal SNPs chr${chr}
set +x
mkdir -p $tmpDir
# calculate recal: VariantRecalibrator
for maxGauLoc in \$(seq 3 6 | sort -r)
do
if [[ -s ${output}_chr${chr}_SNPs_filtered.vcf.gz ]]
then
    break
fi
$java -Djava.io.tmpdir=${tmpDir} -Xmx${memoSmall}g -jar ${GATK} -R $fasta  -L $chr \
   -T VariantRecalibrator \
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
$java -Xmx${memoSmall}g -jar ${GATK} -R $fasta -L $chr \
   -T ApplyRecalibration \
   --ts_filter_level 99.5 \
   --recal_file ${output}_chr${chr}_SNPs_combrec.recal \
   --tranches_file ${output}_chr${chr}_SNPs_combtranch --mode SNP \
   --input ${output}_chr${chr}_SNPs.vcf.gz --out ${output}_chr${chr}_SNPs_filtered.vcf.gz
done " > ${scripts_folder}/subscript_chr${chr}.sh
       	fi
    done
}

##################################################
function recal_indels() {
for chr in `seq 1 22` X
    do
        #### creates the tmpDir if needed
        tmpDir=/scratch0/GATK_chr${chr}
	rm -f ${scripts_folder}/subscript_chr${chr}.sh
        f1=${output}_chr${chr}_indels_hard_filtered.vcf.gz
	f2=${output}_chr${chr}_indels_filtered.vcf.gz
	if [[ ! -s $f1 || ! -s $f2 ]]
        then
        echo "
############### recal_indels chr${chr}
set +x
mkdir -p $tmpDir
#### apply the hard filtering for the indels
$java -Djava.io.tmpdir=${tmpDir} -Xmx${memoSmall}g -jar ${GATK} -R $fasta  -L $chr \
   -T VariantFiltration \
   --filterExpression \"QD < 2.0 || FS > 50.0 || ReadPosRankSum < -20.0\" \
   --filterName \"FAIL\" \
   -V ${output}_chr${chr}_indels.vcf.gz --out ${output}_chr${chr}_indels_hard_filtered.vcf.gz
# calculate recal: VariantRecalibrator
for maxGaussiansIndels in \$(seq 3 6 | sort -r)
do
if [[ -s ${output}_chr${chr}_indels_filtered.vcf.gz ]]
then
    break
fi
$java -Djava.io.tmpdir=${tmpDir} -Xmx${memoSmall}g -jar ${GATK} -R $fasta  -L $chr \
    -T VariantRecalibrator \
   -resource:mills,known=true,training=true,truth=true,prior=12.0 ${bundle}/Mills_and_1000G_gold_standard.indels.b37.vcf \
   -an QD -an FS -an ReadPosRankSum -an InbreedingCoeff \
   -tranche 100.0 -tranche 99.5  -tranche 99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \
   --minNumBadVariants ${numBadIndels} \
   --maxGaussians \${maxGaussiansIndels} \
   --input ${output}_chr${chr}_indels.vcf.gz  --mode INDEL \
   -recalFile ${output}_chr${chr}_indels_combrec.recal \
   -tranchesFile ${output}_chr${chr}_indels_combtranch \
   -rscriptFile ${output}_chr${chr}_recal_plots_indels.R
/share/apps/R/bin/Rscript ${output}_chr${chr}_recal_plots_indels.R
#apply_recal: ApplyRecalibration
$java -Xmx${memoSmall}g -jar ${GATK} -R $fasta  -L $chr \
    -T ApplyRecalibration  \
   --ts_filter_level 95.0 \
   --recal_file ${output}_chr${chr}_indels_combrec.recal \
   --tranches_file ${output}_chr${chr}_indels_combtranch --mode INDEL \
   --input ${output}_chr${chr}_indels.vcf.gz --out ${output}_chr${chr}_indels_filtered.vcf.gz
done " > ${scripts_folder}/subscript_chr${chr}.sh
        fi
    done
}


##################################################
function recal_merge() {
for chr in `seq 1 22` X
    do
        #### creates the tmpDir if needed
        tmpDir=/scratch0/GATK_chr${chr}
        rm -f ${scripts_folder}/subscript_chr${chr}.sh
	f=${output}_chr${chr}_filtered.vcf.gz
        if [[ ! -s $f ]]
        then
        echo "
#### Now we merge SNPs and indels
$java -Djava.io.tmpdir=${tmpDir} -Xmx${memoSmall}g -jar ${GATK} -R $fasta -L $chr \
   -T CombineVariants --assumeIdenticalSamples \
   --variant:SNPs ${output}_chr${chr}_SNPs_filtered.vcf.gz \
   --variant:indels ${output}_chr${chr}_indels_filtered.vcf.gz \
   -genotypeMergeOptions PRIORITIZE  \
   -priority SNPs,indels \
   --out ${output}_chr${chr}_filtered.vcf.gz" > ${scripts_folder}/subscript_chr${chr}.sh
	fi
    done
}

##################################################
function recal() {
    for chr in `seq 1 22` X
    do
        f=${output}_chr${chr}_filtered.vcf
        if [[ ! -s $f || "$force" == "yes" ]]
        then 
            #### creates the tmpDir if needed
            tmpDir=/scratch0/GATK_chr${chr}
            maxGauLoc=${maxGaussians}
            if [[ "$chr" == "14" ]]
            then
                maxGauLoc=4
                maxGaussiansIndels=4
            fi
            echo "
############### recal chr${chr}
set +x
mkdir -p $tmpDir
####### first SNPs
#### extract the SNPs: SelectVariants 
$java  -Djava.io.tmpdir=${tmpDir} -Xmx${memoSmall}g -jar ${GATK} -R $fasta -L $chr \
   -T SelectVariants \
   -selectType SNP \
   -V ${output}_chr${chr}.vcf.gz  --out ${output}_chr${chr}_SNPs.vcf.gz
# calculate recal: VariantRecalibrator 
$java -Djava.io.tmpdir=${tmpDir} -Xmx${memoSmall}g -jar ${GATK} -R $fasta  -L $chr \
   -T VariantRecalibrator \
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
$java -Xmx${memoSmall}g -jar ${GATK} -R $fasta -L $chr \
   -T ApplyRecalibration \
   --ts_filter_level 99.5 \
   --recal_file ${output}_chr${chr}_SNPs_combrec.recal \
   --tranches_file ${output}_chr${chr}_SNPs_combtranch --mode SNP \
   --input ${output}_chr${chr}_SNPs.vcf.gz --out ${output}_chr${chr}_SNPs_filtered.vcf.gz 
####### now indels
#### extract the indels
$java  -Djava.io.tmpdir=${tmpDir} -Xmx${memoSmall}g -jar ${GATK}  -R $fasta  -L $chr \
   -T SelectVariants \
   -selectType INDEL \
   -selectType MIXED \
   -V ${output}_chr${chr}.vcf.gz --out ${output}_chr${chr}_indels.vcf.gz
#### apply the hard filtering for the indels
$java -Djava.io.tmpdir=${tmpDir} -Xmx${memoSmall}g -jar ${GATK} -R $fasta  -L $chr \
   -T VariantFiltration \
   --filterExpression \"QD < 2.0 || FS > 50.0 || ReadPosRankSum < -20.0\" \
   --filterName \"FAIL\" \
   -V ${output}_chr${chr}_indels.vcf.gz --out ${output}_chr${chr}_indels_hard_filtered.vcf.gz
# calculate recal: VariantRecalibrator 
$java -Djava.io.tmpdir=${tmpDir} -Xmx${memoSmall}g -jar ${GATK} -R $fasta  -L $chr \
    -T VariantRecalibrator \
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
$java -Xmx${memoSmall}g -jar ${GATK} -R $fasta  -L $chr \
    -T ApplyRecalibration  \
   --ts_filter_level 95.0 \
   --recal_file ${output}_chr${chr}_indels_combrec.recal \
   --tranches_file ${output}_chr${chr}_indels_combtranch --mode INDEL \
   --input ${output}_chr${chr}_indels.vcf.gz --out ${output}_chr${chr}_indels_filtered.vcf.gz
#### Now we merge SNPs and indels
$java -Djava.io.tmpdir=${tmpDir} -Xmx${memoSmall}g -jar ${GATK} -R $fasta -L $chr \
   -T CombineVariants --assumeIdenticalSamples \
   --variant:SNPs ${output}_chr${chr}_SNPs_filtered.vcf.gz \
   --variant:indels ${output}_chr${chr}_indels_filtered.vcf.gz \
   -genotypeMergeOptions PRIORITIZE  \
   -priority SNPs,indels \
   --out ${output}_chr${chr}_filtered.vcf.gz
rm -rf $tmpDir
rm ${output}_chr${chr}_SNPs.vcf.gz*
rm ${output}_chr${chr}_SNPs_filtered.vcf.gz*
rm ${output}_chr${chr}_indels.vcf.gz*
# keep indels for Costin
# rm ${output}_chr${chr}_indels_filtered.vcf.gz*
# rm ${output}_chr${chr}_indels_hard_filtered.vcf.gz*
" >> ${scripts_folder}/subscript_chr${chr}.sh
    else
        echo $f exists
	fi
    done
}

function file_exists() {
    if [[ ! -s $1 ]]
    then
        error "$1 does not exist"
    fi
}

function check_recal() {
    for chr in `seq 1 22` X
    do
        file_exists ${output}_chr${chr}_filtered.vcf.gz
        file_exists ${output}_chr${chr}_filtered.vcf.gz.tbi
        file_exists ${output}_chr${chr}_indels_combrec.recal
        file_exists ${output}_chr${chr}_indels_combrec.recal.idx
        file_exists ${output}_chr${chr}_indels_combtranch
        file_exists ${output}_chr${chr}_indels_filtered.vcf.gz
        file_exists ${output}_chr${chr}_indels_filtered.vcf.gz.tbi
        file_exists ${output}_chr${chr}_indels_hard_filtered.vcf.gz
        file_exists ${output}_chr${chr}_indels_hard_filtered.vcf.gz.tbi
        file_exists ${output}_chr${chr}_recal_plots_indels.R
        file_exists ${output}_chr${chr}_recal_plots_indels.R.pdf
        file_exists ${output}_chr${chr}_recal_plots_snps.R
        file_exists ${output}_chr${chr}_recal_plots_snps.R.pdf
        file_exists ${output}_chr${chr}_SNPs_combrec.recal
        file_exists ${output}_chr${chr}_SNPs_combrec.recal.idx
        file_exists ${output}_chr${chr}_SNPs_combtranch
        file_exists ${output}_chr${chr}.vcf.gz
        file_exists ${output}_chr${chr}.vcf.gz.tbi
    done
}
    


##################################################
function annovar() {
    check_recal
    memo=15
    for chr in `seq 1 22` X
    do
        if [[ ! -s ${output}_snpStats/annotations_chr${chr}.csv  || "$force" == "yes" ]]
        then
            echo "
############### annovar chr${chr}
# creates ${output}_chr${chr}_recal_filtered2.vcf 
zcat ${output}_chr${chr}_filtered.vcf.gz | perl ${baseFolder}/UCLex/custom_filtering.pl ${output}_chr${chr}_filtered2.vcf ${GQ}
cut -f1-8 ${output}_chr${chr}_filtered2.vcf > ${output}_chr${chr}_for_annovar.vcf
/cluster/project8/vyp/vincent/Software_heavy/annovar_Feb2013/convert2annovar.pl --allallele -format vcf4 --includeinfo ${output}_chr${chr}_for_annovar.vcf > ${output}_chr${chr}_db
/cluster/project8/vyp/vincent/Software_heavy/annovar_Feb2013/summarize_annovar_VP.pl -ver1000g 1000g2012apr -verdbsnp 137 -veresp 6500si -alltranscript -buildver hg19 --genetype ensgene --remove ${output}_chr${chr}_db /cluster/project8/vyp/vincent/Software_heavy/annovar_Feb2013/humandb_hg19/ 
python ${baseFolder}/UCLex/annovar_vcf_combine_VP.py ${output}_chr${chr}_filtered2.vcf ${output}_chr${chr}_db.exome_summary.csv ${output}_chr${chr}_exome_table.csv
perl ${baseFolder}/UCLex/make_matrix_calls.pl ${output}_chr${chr}_exome_table.csv ${output} $chr
# compress ${output}_chr${chr}_filtered2.vcf 
/share/apps/genomics/htslib-1.1/bin/bgzip -f -c ${output}_chr${chr}_filtered2.vcf > ${output}_chr${chr}_filtered2.vcf.gz
/share/apps/genomics/htslib-1.1/bin/tabix -f -p vcf ${output}_chr${chr}_filtered2.vcf.gz
rm ${output}_chr${chr}_filtered2.vcf
# compress ${output}_chr${chr}_for_annovar.vcf
/share/apps/genomics/htslib-1.1/bin/bgzip -f -c ${output}_chr${chr}_for_annovar.vcf > ${output}_chr${chr}_for_annovar.vcf.gz
/share/apps/genomics/htslib-1.1/bin/tabix -f -p vcf ${output}_chr${chr}_for_annovar.vcf.gz
rm ${output}_chr${chr}_for_annovar.vcf
" >> ${scripts_folder}/subscript_chr${chr}.sh
    else
        echo ${output}_snpStats/annotations_chr${chr}.csv exists
    fi
    done
}



##################################################
# write straight to mongo
function VEP_mongo() {
    memo=30
    mkdir -p ${output}_VEP
    #VEP_output="--tab --output_file ${output}_VEP/VEP_${chr}.txt"
    ensembl=/cluster/project8/vyp/AdamLevine/software/ensembl/
    VEP_DIR=/cluster/project8/vyp/Software/ensembl-tools-release-82/scripts/
    DIR_CACHE=/SAN/vyplab/NCMD_raw/VEP/cache/
    DIR_PLUGINS=${DIR_CACHE}/Plugins
    for chr in `seq 1 22` X
    do
        #output_lines=`zcat ${output}_VEP/VEP_${chr}.json.gz | wc -l` 
        #input_lines=`tail -n+2 ${output}_VEP/chr${chr}_for_VEP.vcf | wc -l`
        #echo ${output}_VEP/chr${chr}_for_VEP.vcf input lines: $input_lines
        #echo ${output}_VEP/VEP_${chr}.json.gz output lines: $output_lines
        #if [[ $lines -gt 0 ]]; then echo ${output}_VEP/VEP_${chr}.json.gz $lines gt than 0, skipping; continue; fi
        echo "
############### VEP_mongo chr${chr}
# split single lines
zcat ${output}_chr${chr}_for_annovar.vcf.gz | /share/apps/python/bin/python ${baseFolder}/annotation/multiallele_to_single_vcf.py --headers CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO > ${output}_VEP/chr${chr}_for_VEP.vcf
/share/apps/genomics/htslib-1.1/bin/bgzip -f -c ${output}_VEP/chr${chr}_for_VEP.vcf > ${output}_VEP/chr${chr}_for_VEP.vcf.gz
/share/apps/genomics/htslib-1.1/bin/tabix -f -p vcf ${output}_VEP/chr${chr}_for_VEP.vcf.gz
rm ${output}_VEP/chr${chr}_for_VEP.vcf
#### RUN CADD
/share/apps/genomics/CADD_v1.3/bin/score.sh ${output}_VEP/chr${chr}_for_VEP.vcf.gz ${output}_VEP/CADD_chr${chr}.vcf.gz
/share/apps/genomics/htslib-1.1/bin/tabix -p vcf ${output}_VEP/CADD_chr${chr}.vcf.gz
####CONFIGURE SOFTWARE SHORTCUTS AND PATHS
reference=1kg
ensembl=/cluster/project8/vyp/AdamLevine/software/ensembl/
export PERL5LIB=${PERL5LIB}:${ensembl}/src/bioperl-1.6.1::${ensembl}/src/ensembl/modules:${ensembl}/src/ensembl-compara/modules:${ensembl}/src/ensembl-variation/modules:${ensembl}/src/ensembl-funcgen/modules:${ensembl}/Plugins
export PATH=$PATH:/cluster/project8/vyp/vincent/Software/tabix-0.2.5/
# RUN VEP
/share/apps/perl/bin/perl ${VEP_DIR}/variant_effect_predictor/variant_effect_predictor.pl --offline \
--ASSEMBLY GRCh37 --fasta /SAN/vyplab/UKIRDC/reference/human_g1k_v37.fasta \
--cache --dir_cache ${DIR_CACHE} \
--dir_plugins ${DIR_PLUGINS} \
--sift b --polyphen b --symbol --canonical --check_existing --check_alleles  \
--fork 4 --maf_esp --gmaf --maf_1kg --maf_exac \
--no_progress --quiet \
--custom /cluster/project9/IBDAJE/VEP_custom_annotations/1kg/1kg/chr${chr}_EUR.vcf.gz,1KG_EUR,vcf,exact \
--custom /cluster/project9/IBDAJE/VEP_custom_annotations/1kg/1kg/chr${chr}_AFR.vcf.gz,1KG_AFR,vcf,exact \
--custom /cluster/project9/IBDAJE/VEP_custom_annotations/1kg/1kg/chr${chr}_AMR.vcf.gz,1KG_AMR,vcf,exact \
--custom /cluster/project9/IBDAJE/VEP_custom_annotations/1kg/1kg/chr${chr}_ASN.vcf.gz,1KG_ASN,vcf,exact \
--custom /cluster/project9/IBDAJE/VEP_custom_annotations/1kg/esp/chr${chr}_EA.vcf.gz,ESP_EA,vcf,exact \
--custom /cluster/project9/IBDAJE/VEP_custom_annotations/1kg/esp/chr${chr}_AA.vcf.gz,ESP_AA,vcf,exact \
--custom /cluster/scratch3/vyp-scratch2/reference_datasets/Kaviar/Kaviar-160204-Public/hg19/VEP_annotation.vcf.gz,Kaviar,vcf,exact \
--custom ${output}_VEP/CADD_chr${chr}.vcf.gz,CADD,vcf,exact \
--plugin Condel,${DIR_PLUGINS}/config/Condel/config,b \
--plugin Carol \
--no_stats \
--hgvs \
--pubmed \
--plugin HGVSshift \
--plugin SameCodon \
--input_file ${output}_VEP/chr${chr}_for_VEP.vcf.gz \
--json \
--output_file STDOUT | python ${baseFolder}/annotation/postprocess_VEP_json.py --release mainset_${currentUCLex} | grep '^JSON:' | sed 's/^JSON://' > ${output}_VEP/VEP_chr${chr}.json
" >> ${scripts_folder}/subscript_chr${chr}.sh
#--output_file STDOUT | python ${baseFolder}/annotation/postprocess_VEP_json.py | grep '^JSON:' | sed 's/^JSON://' | mongoimport --db uclex --collection variants --host phenotips
  done
}


##################################################
# 
function VEP() {
    memo=30
    mkdir -p ${output}_VEP
    #VEP_output="--tab --output_file ${output}_VEP/VEP_${chr}.txt"
    ensembl=/cluster/project8/vyp/AdamLevine/software/ensembl/
    VEP_DIR=/cluster/project8/vyp/Software/ensembl-tools-release-82/scripts/
    DIR_CACHE=/SAN/vyplab/NCMD_raw/VEP/cache/
    DIR_PLUGINS=${DIR_CACHE}/Plugins
    for chr in `seq 1 22` X
    do
        #output_lines=`zcat ${output}_VEP/VEP_${chr}.json.gz | wc -l` 
        #input_lines=`tail -n+2 ${output}_VEP/chr${chr}_for_VEP.vcf | wc -l`
        #echo ${output}_VEP/chr${chr}_for_VEP.vcf input lines: $input_lines
        #echo ${output}_VEP/VEP_${chr}.json.gz output lines: $output_lines
        #if [[ $lines -gt 0 ]]; then echo ${output}_VEP/VEP_${chr}.json.gz $lines gt than 0, skipping; continue; fi
        echo "
############### VEP chr${chr}
# split single lines
zcat ${output}_chr${chr}_filtered.vcf.gz | /share/apps/python/bin/python ${baseFolder}/annotation/multiallele_to_single_vcf.py | gzip > ${output}_chr${chr}_filtered3.vcf.gz 
# make genotype matrix
#mainset_July2016_chr${chr}_filtered3-genotypes.csv
#mainset_July2016_chr${chr}_filtered3-annotations.csv
#mainset_July2016_chr${chr}_filtered3-genotypes_depth.csv
/share/apps/python/bin/python ${baseFolder}/annotation/postprocess_VEP_vcf.py --file ${output}_chr${chr}_filtered3.vcf.gz --genotypes --depth
rm ${output}_chr${chr}_filtered3-annotations.csv
gzip -f ${output}_chr${chr}_filtered3-genotypes.csv
gzip -f ${output}_chr${chr}_filtered3-genotypes_depth.csv
zcat ${output}_chr${chr}_filtered3.vcf.gz | /share/apps/python/bin/python ${baseFolder}/annotation/multiallele_to_single_vcf.py --headers CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT | cut -f1-9 | gzip > ${output}_VEP/chr${chr}_for_VEP.vcf.gz
#### RUN CADD
/share/apps/genomics/CADD_v1.3/bin/score.sh ${output}_VEP/chr${chr}_for_VEP.vcf.gz ${output}_VEP/CADD_chr${chr}.vcf.gz
/share/apps/genomics/htslib-1.1/bin/tabix -p vcf ${output}_VEP/CADD_chr${chr}.vcf.gz
####CONFIGURE SOFTWARE SHORTCUTS AND PATHS
reference=1kg
export PERL5LIB=${PERL5LIB}:${ensembl}/src/bioperl-1.6.1::${ensembl}/src/ensembl/modules:${ensembl}/src/ensembl-compara/modules:${ensembl}/src/ensembl-variation/modules:${ensembl}/src/ensembl-funcgen/modules:${ensembl}/Plugins
export PATH=$PATH:/cluster/project8/vyp/vincent/Software/tabix-0.2.5/
# RUN VEP
/share/apps/perl/bin/perl ${VEP_DIR}/variant_effect_predictor/variant_effect_predictor.pl --offline \
--ASSEMBLY GRCh37 --fasta /SAN/vyplab/UKIRDC/reference/human_g1k_v37.fasta \
--cache --dir_cache ${DIR_CACHE} \
--dir_plugins ${DIR_PLUGINS} \
--sift b --polyphen b --symbol --canonical --check_existing --check_alleles  \
--fork 4 --maf_esp --gmaf --maf_1kg --maf_exac \
--pubmed \
--no_progress --quiet \
--custom /cluster/scratch3/vyp-scratch2/reference_datasets/Kaviar/Kaviar-160204-Public/hg19/VEP_annotation.vcf.gz,Kaviar,vcf,exact \
--custom ${output}_VEP/CADD_chr${chr}.vcf.gz,CADD,vcf,exact \
--plugin Condel,${DIR_PLUGINS}/config/Condel/config,b \
--plugin Carol \
--no_stats \
--hgvs \
--plugin HGVSshift \
--plugin SameCodon \
--input_file ${output}_VEP/chr${chr}_for_VEP.vcf.gz \
--pick \
--everything \
--tab \
--fields Uploaded_variation,Location,Allele,Gene,Feature,Feature_type,Consequence,cDNA_position,CDS_position,Protein_position,Amino_acids,Codons,Existing_variation,IMPACT,DISTANCE,STRAND,FLAGSVARIANT_CLASS,SYMBOL,SYMBOL_SOURCE,HGNC_ID,BIOTYPE,CANONICAL,TSL,APPRIS,CCDS,ENSP,SWISSPROT,TREMBL,UNIPARC,GENE_PHENO,SIFT,PolyPhen,EXON,INTRON,DOMAINS,HGVSc,HGVSp,HGVS_OFFSET,GMAF,AFR_MAF,AMR_MAF,EAS_MAF,EUR_MAF,SAS_MAF,AA_MAF,EA_MAF,ExAC_MAF,ExAC_Adj_MAF,ExAC_AFR_MAF,ExAC_AMR_MAF,ExAC_EAS_MAF,ExAC_FIN_MAF,ExAC_NFE_MAF,ExAC_OTH_MAF,ExAC_SAS_MAF,CLIN_SIG,SOMATIC,PHENO,PUBMED,MOTIF_NAME,MOTIF_POS,HIGH_INF_POS,MOTIF_SCORE_CHANGE,CAROL,HGVSc_unshifted,HGVSp_unshifted,SameCodon,Kaviar \
--output_file STDOUT > ${output}_VEP/VEP_chr${chr}.tab
/share/apps/R/bin/Rscript ${baseFolder}/annotation/postprocess_VEP_tab.R --input ${output}_VEP/VEP_chr${chr}.tab --output ${output}_VEP/VEP_chr${chr}.csv
gzip -f ${output}_VEP/VEP_chr${chr}.csv
gzip -f ${output}_VEP/VEP_chr${chr}.tab
" >> ${scripts_folder}/subscript_chr${chr}.sh
  done
}


function check_annovar() {
    for chr in `seq 1 22` X
    do
        file_exists ${output}_snpStats/annotations_chr${chr}.csv
    done
}

##################################################
function convertToR() {
    check_annovar
    mkdir -p ${output}_snpStats
    for chr in `seq 1 22` X
    do
        if [[ ! -s ${output}_snpStats/chr${chr}_snpStats.RData || "$force" == "yes" ]]
        then
            echo "
############### convertToR chr${chr}
/share/apps/R/bin/Rscript ${baseFolder}/UCLex/convert_to_R.R --chromosome ${chr} --root ${output} > ${output}_snpStats/convert_to_R_chr${chr}.out
" >> ${scripts_folder}/subscript_chr${chr}.sh
        fi
    done
}

##################################################
function finalCrunch() {
    keyWords=data/controlKeywords.tab
    casekeyWords=data/caseKeywords.tab
    for chr in `seq 1 22` X
    do
        script=cluster/submission/subscript_chr${chr}.sh
        echo "
${baseFolder}/UCLex/crunch_controls.pl ${output}_chr${chr}_exome_table.csv $keyWords $casekeyWords ${output}_chr${chr}_exome_crunched.csv data/sampleList_exome.tab none no  ##include all samples
" >> $script
    done
}


##################################################
function all() {
    for chr in `seq 1 22` X
    do
        script=${scripts_folder}/subscript_chr${chr}.sh
        echo "
#! /bin/bash
set +x
        " > ${script}
        chmod u+rx ${script}
    done
    genotype
    recal
    annovar
    convertToR
    h_vmem=20
    t_mem=20
}


for m in `echo ${mode} | tr ',' '\n'`
do
    echo mode ${m}
    ${m}
done

job_n=$(ls ${scripts_folder}/subscript_chr* | wc -l)
if [[ ${job_n} -gt 0 ]]
then
    echo "
#$ -o ${stdout_folder}
#$ -e ${stderr_folder}
#$ -S /bin/bash
#$ -l h_vmem=${memo}G,tmem=${memo}G
#$ -l h_rt=240:0:0
#$ -cwd 
#$ -t 1-${job_n}

LISTJOBS=(dud $(ls ${scripts_folder}/subscript_chr* | sort -V | tr '\n' ' '))

CHR=\${LISTJOBS[ \$SGE_TASK_ID ]}

echo \$HOSTNAME
hostname
date

sh \${CHR}

date

" > $mainScript
else
    rm ${scripts_folder}/calling.sh
    echo "No jobs to submit!"
fi



##############################
for chr in `seq 1 22` X
do
    ls -ltrh ${scripts_folder}/subscript_chr${chr}.sh
done

