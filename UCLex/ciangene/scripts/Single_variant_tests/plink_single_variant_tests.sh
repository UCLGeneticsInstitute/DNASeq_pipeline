##!/bin/bash
shopt -s expand_aliases
source ~/.bashrc

R=/share/apps/R-3.1.0/bin/R
runSh='sh /cluster/project8/vyp/cian/scripts/bash/runBashCluster_large.sh'
rootODir=$1
release=$2
#rootODir=/scratch2/vyp-scratch2/ciangene
rootODir=/scratch2/vyp-scratch2/cian
release=July2015
rootODir=${1-$rootODir}
release=${2-$release}
bDir=${rootODir}/UCLex_${release}/

plink=/home/sejjcmu/bin/plink/plink
Groups=${bDir}cohort.list
nbGroups=$(wc -l  $Groups | awk '{print $1}') 
data=$bDir"allChr_snpStats_out"
Pheno=${bDir}Clean_pheno_subset_plink
permPheno=${bDir}Phenotypes_fastlmm_permuted
covar=$bDir"TechPCs.vect"
phenos=$(wc -l  $Groups | awk '{print $1}') 

oDir=$bDir"Single_variant_tests/"
if [ ! -e $oDir ]; then mkdir $oDir; fi

cwd=$(pwd)
for pheno in $(seq 1 16)
do
	batch=$(sed -n $pheno'p' $Groups)	
	oFile=$oDir"/run_${batch}.sh"
	echo "
	$plink --noweb  --bfile $data --assoc --counts --allow-no-sex  --pheno $Pheno \
	--adjust --mpheno $pheno --out $oDir$batch"_counts_assoc"
	#$plink --noweb  --bfile $data --fisher --allow-no-sex  --pheno $Pheno --adjust --mpheno $pheno --out $oDir$batch"_fisher"
	#$plink --noweb  --bfile $data --logistic --allow-no-sex --pheno $Pheno --adjust --mpheno $pheno --out $oDir$batch"_logistic_tk_depth"
   # $plink --noweb  --bfile $data --fisher --allow-no-sex  --pheno $permPheno --adjust --mpheno $pheno --out $oDir$batch"_fisher_perm"	
	##$plink --noweb  --bfile $data --logistic  --hide-covar --allow-no-sex --pheno $Pheno --adjust --mpheno $pheno --covar $covar --covar-number 1-2 --out $oDir$batch"_logistic_tech_pcs_covars"
	" > $oFile
	cd $oDir; $runSh "run_${batch}.sh" ; cd $cwd
done






