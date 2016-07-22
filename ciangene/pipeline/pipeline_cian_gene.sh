Rbin=/cluster/project8/vyp/vincent/Software/R-3.1.2/bin/R

## First step - genotype matrix and initial filtering
firstStep=${repo}/scripts/first.step.R  ##step 1.1
clean=${repo}/scripts/variant_filtering/qc.sh   ##step 1.2
filter=${repo}/scripts/variant_filtering/filter_snps.R   ##step 1.3
pca=${repo}/scripts/PCA/ancestry_pca.R  ## step 1.4
pca_extract=${repo}/scripts/PCA/getPCAsnps_UCLex.sh
plot_pca=${repo}/scripts/PCA/plot_pca.R


## second step
pheno=${repo}/scripts/MakePhenotypes/prepare_all_phenos.R ## step 2.1
MissingNess=${repo}/scripts/CaseControlMissingness.R ##step 2.2

## Third step - creating and validating the technical Kinship
VariantLists=${repo}/scripts/LDAK/make_list_of_variants_for_gene_tests.R
secondStep=${repo}/scripts/convert_genotype_to_missingNonMissing.sh  ## step2
makeKin=${repo}/scripts/make_kinships_new.sh # step 2.1
checkKin=${repo}/scripts/check_tk_residuals.sh # step 2.2
convertKin=${repo}/scripts/prepare_kinship_matrix_for_fastLMM.R # step 2.3

## Fourth step - case control tests
thirdStep=${repo}/scripts/LDAK/run_ldak_on_all_phenos.sh
singleVariant=${repo}/scripts/plink_single_variant_tests.sh
fastlmm=${repo}/scripts/run_fastLMM_on_all_phenotypes.sh

##### default value of all arguments
rootODir=/scratch2/vyp-scratch2/cian/
release=February2015

######## create directories
for dir in cluster cluster/submission cluster/error cluster/out cluster/R; do
    if [ ! -e $dir ]; then mkdir $dir; fi
done




until [ -z "$1" ]; do
    # use a case statement to test vars. we always test $1 and shift at the end of the for block.
    case $1 in
	--rootODir)
	    shift
	    rootODir=$1;;
	--script)
	    shift
	    script=$1;;
	--release)
	    shift
	    release=$1;;
	--step1)
	    shift
	    step1=$1;;
	--step2)
	    shift
	    step2=$1;;
	--step3)
	    shift
	    step3=$1;;
	--step4)
	    shift
	    step4=$1;;
	-* )
	    echo "Unrecognized option: $1"
	    exit 1;;
    esac
    shift
    if [ "$#" = "0" ]; then break; fi
    echo $1
done


hold=""


##########
if [[ "$step1" == "yes" ]]; then    

    script=cluster/submission/step1.sh

    echo "
#$ -S /bin/bash
#$ -o cluster/out
#$ -e cluster/error
#$ -l h_vmem=10G,tmem=10G
#$ -pe smp 1
#$ -N step1_cian
#$ -l h_rt=24:00:00
#$ -cwd

$Rbin CMD BATCH --no-save --no-restore --release=${release} --rootODir=${rootODir} $firstStep cluster/R/step1.1_first_step.Rout
sh $clean $rootODir $release
$Rbin CMD BATCH --no-save --no-restore --release=${release} --rootODir=${rootODir} $filter cluster/R/step1.3.filter_snps.Rout
$Rbin CMD BATCH --no-save --no-restore --release=${release} --rootODir=${rootODir} $pca cluster/R/step1.4.pca.Rout
sh $pca_extract $rootODir $release
$Rbin CMD BATCH --no-save --no-restore --release=${release} --rootODir=${rootODir} $plot_pca cluster/R/step1.4.Plotpca.Rout


" > $script
    
    qsub $hold $script
    if [[ "$hold" == "" ]]; then hold="-hold_jid step1_cian"; else hold="$hold,step1_cian"; fi
fi



#########################################
if [[ "$step2" == "yes" ]]; then    

    script=cluster/submission/step2.sh

    echo "
#$ -S /bin/bash
#$ -o cluster/out
#$ -e cluster/error
#$ -l h_vmem=2G,tmem=2G
#$ -pe smp 4
#$ -N step2_cian
#$ -l h_rt=24:00:00
#$ -cwd

$Rbin CMD BATCH --no-save --no-restore --release=${release} --rootODir=${rootODir} $pheno cluster/R/step2.1.pheno.Rout
$Rbin CMD BATCH --no-save --no-restore --release=${release} --rootODir=${rootODir} $MissingNess cluster/R/step2.2.CaseControlMissingness.Rout

" > $script
    
    qsub $hold $script
    if [[ "$hold" == "" ]]; then hold="-hold_jid step2_cian"; else hold="$hold,step2_cian"; fi
fi

########################################

#########
if [[ "$step3" == "yes" ]]; then    

    script=cluster/submission/step3.sh

    echo "
#$ -S /bin/bash
#$ -o cluster/out
#$ -e cluster/error
#$ -l h_vmem=5G,tmem=5G
#$ -pe smp 1
#$ -N step3_cian
#$ -l h_rt=24:00:00
#$ -cwd

$Rbin CMD BATCH --no-save --no-restore --release=${release} --rootODir=${rootODir} $VariantLists cluster/R/step3.1.variant.lists.Rout
sh $secondStep $rootODir $release ## convert geno to missingNonMissing
sh $makeKin $rootODir $release ### make kinship matrices
sh $checkKin $rootODir $release # check how much variance the kinships explained. 
$Rbin CMD BATCH --no-save --no-restore --release=${release} --rootODir=${rootODir} $convertKin cluster/R/step2.Rout
" > $script

    qsub $hold $script
    if [[ "$hold" == "" ]]; then hold="-hold_jid step3_cian"; else hold="$hold,step3_cian"; fi

fi



#########
if [[ "$step4" == "yes" ]]; then    

    script=cluster/submission/step4.sh

    echo "
#$ -S /bin/bash
#$ -o cluster/out
#$ -e cluster/error
#$ -l h_vmem=40G,tmem=40G
#$ -pe smp 1
#$ -N step4_cian
#$ -l h_rt=24:00:00
#$ -cwd

sh $thirdStep $rootODir $release 

sh $singleVariant $rootODir $release 

sh $fastlmm $rootODir $release 

" > $script

    qsub $hold $script
    if [[ "$hold" == "" ]]; then hold="-hold_jid step4_cian"; else hold="$hold,step4_cian"; fi

fi



ls -ltrh $script

