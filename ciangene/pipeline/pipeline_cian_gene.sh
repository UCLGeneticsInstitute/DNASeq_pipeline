Rbin=/share/apps/R/bin/R

## First step - genotype matrix and initial filtering
clean=${repo}/variant_filtering/qc.sh   ##step 1.2
filter=${repo}/variant_filtering/filter_snps.R   ##step 1.3
pca=${repo}/PCA/ancestry_pca.R  ## step 1.4
pca_extract=${repo}/PCA/getPCAsnps_UCLex.sh
plot_pca=${repo}/PCA/plot_pca.R


## Third step - creating and validating the technical Kinship
VariantLists=${repo}/Gene_Based_Tests/make_list_of_variants_for_gene_tests.R # step 3.1
secondStep=${repo}/convert_genotype_to_missingNonMissing.sh  ## step3.2
makeKin=${repo}/Make_Kinships/make_kinships_new.sh # step 3.3
checkKin=${repo}/check_tk_residuals.sh # step 3.4
convertKin=${repo}/prepare_kinship_matrix_for_fastLMM.R # step 3.5

## Fourth step - case control tests
thirdStep=${repo}/LDAK/run_ldak_on_all_phenos.sh
singleVariant=${repo}/plink_single_variant_tests.sh
fastlmm=${repo}/run_fastLMM_on_all_phenotypes.sh


## Fifth step - ExomeDepth


##### default value of all arguments
#rootODir=/cluster/scratch3/vyp-scratch2/cian
rootODir=/cluster/project8/vyp/cian/data/UCLex
release=June2016

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
	    step4=$1;;	--step5)
	    shift
	    step5=$1;;
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
#$ -l h_vmem=20G,tmem=20G
#$ -pe smp 1
#$ -N step1_cian
#$ -l h_rt=24:00:00
#$ -cwd

$Rbin CMD BATCH --no-save --no-restore --release=${release} --rootODir=${rootODir} ${repo}/first.step.R cluster/R/step1.1_first_step.Rout
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
#$ -l h_vmem=25G,tmem=25G
#$ -pe smp 1
#$ -N step2_cian
#$ -l h_rt=24:00:00
#$ -cwd

$Rbin CMD BATCH --no-save --no-restore --release=${release} --rootODir=${rootODir} ${repo}/checkSteps.R                             cluster/R/checkSteps.Rout
$Rbin CMD BATCH --no-save --no-restore --release=${release} --rootODir=${rootODir} ${repo}/MakePhenotypes/make_phenotype_file.R     cluster/R/make_phenotype_file.Rout
$Rbin CMD BATCH --no-save --no-restore --release=${release} --rootODir=${rootODir} ${repo}/MakePhenotypes/CaseControl_support.R     cluster/R/CaseControl_support.Rout
$Rbin CMD BATCH --no-save --no-restore --release=${release} --rootODir=${rootODir} ${repo}/MakePhenotypes/prepare_all_phenos.R      cluster/R/step2.1.pheno.Rout
$Rbin CMD BATCH --no-save --no-restore --release=${release} --rootODir=${rootODir} ${repo}/MakePhenotypes/MakeGoodPhenotypeFile.R   cluster/R/MakeGoodPhenotypeFile.Rout
$Rbin CMD BATCH --no-save --no-restore --release=${release} --rootODir=${rootODir} ${repo}/CaseControlMissingness.R                 cluster/R/step2.2.CaseControlMissingness.Rout
$Rbin CMD BATCH --no-save --no-restore --release=${release} --rootODir=${rootODir} ${repo}/getReadDepthByGene.R                     cluster/R/getReadDepthByGene.Rout


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
$Rbin CMD BATCH --no-save --no-restore --release=${release} --rootODir=${rootODir} $convertKin cluster/R/step3.2.Rout
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
#$ -l h_vmem=10G,tmem=10G
#$ -pe smp 1
#$ -N step4_cian
#$ -l h_rt=24:00:00
#$ -cwd

#$Rbin CMD BATCH --no-save --no-restore --release=${release} --rootODir=${rootODir} ${repo}/Gene_Based_Tests/make_variant_list_gene_tests.R cluster/R/make_variant_list_gene_tests.Rout
sh ${repo}/Gene_Based_Tests/extract_variants.sh $rootODir $release 
#$Rbin CMD BATCH --no-save --no-restore --release=${release} --rootODir=${rootODir} ${repo}/Gene_Based_Tests/SKAT_uclex.R cluster/R/SKAT_uclex.Rout

#sh $thirdStep $rootODir $release 
#sh $singleVariant $rootODir $release 
#sh $fastlmm $rootODir $release 

" > $script

    qsub $hold $script
    if [[ "$hold" == "" ]]; then hold="-hold_jid step4_cian"; else hold="$hold,step4_cian"; fi

fi


#########
if [[ "$step5" == "yes" ]]; then    

    script=cluster/submission/step5.sh

    echo "
#$ -S /bin/bash
#$ -o cluster/out
#$ -e cluster/error
#$ -l h_vmem=10G,tmem=10G
#$ -pe smp 1
#$ -N step5_cian
#$ -l h_rt=24:00:00
#$ -cwd

#sh ${repo}/ExomeDepth/makeMergedBed.sh $rootODir $release 
$Rbin CMD BATCH --no-save --no-restore --release=${release} --rootODir=${rootODir} ${repo}/ExomeDepth/exomeDepth_rscript.R cluster/R/exomeDepth_rscript.Rout

" > $script

    qsub $hold $script
    if [[ "$hold" == "" ]]; then hold="-hold_jid step5_cian"; else hold="$hold,step5_cian"; fi

fi





ls -ltrh $script
