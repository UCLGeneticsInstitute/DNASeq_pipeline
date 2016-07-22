#/bin/sh

shopt -s expand_aliases
source ~/.bashrc

release=May2015
bDir=/scratch2/vyp-scratch2/cian/UCLex_${release}/
Groups=$bDir"GroupNames"
nbGroups=$(wc -l  $Groups | awk '{print $1}') 

echo $nbGroups "phenotypes for fastLMM to chew on"

oFolder=$bDir"FastLMM_Single_Variant_all_phenos/"
if [ ! -e $oFolder ]; then mkdir $oFolder; fi

runSh='sh /cluster/project8/vyp/cian/scripts/bash/runBashCluster.sh'

for pheno in $(seq 73 $nbGroups)
do

batch=$(sed -n $pheno'p' $Groups)
target=$oFolder"/"$batch".fastlmm.sh"; if [ -e $target ]; then rm $target; fi
oFile=$oFolder$batch"_Out"

echo "
	release=February2015

	bDir=/scratch2/vyp-scratch2/cian/UCLex_${release}/

	fastlmm='/share/apps/bahler/FaSTLMM.207/Bin/Linux_MKL/fastlmmc'
	snpFile=$bDir"allChr_snpStats_out"
	kinshipFile=$bDir"TechnicalKinship_Fastlmm"
	phenoFile=$bDir"Phenotypes_fastlmm"
	extract=$bDir"Clean_variants_Func"
	permPheno=$bDir"Phenotypes_fastlmm_permuted"
	resPheno=$bDir"phenotype_res"
	for chunk in {1..20}
	do
	"'$fastlmm'"  -simLearnType Full -verboseOutput -maxThreads 1 -bfile "'$snpFile'" -sim "'$kinshipFile'" -pheno "'$phenoFile'" -mpheno "$pheno" -out "$batch"_Tech_kin_"'$chunk'" -numjobs 20 -thisjob "'$chunk'"
	"'$fastlmm'"  -linreg -simLearnType Full -verboseOutput -maxThreads 1 -bfile "'$snpFile'" -pheno "'$phenoFile'" -mpheno "$pheno" -out "$batch"_no_kin_"'$chunk'" -numjobs 20 -thisjob "'$chunk'"
	"'$fastlmm'"  -linreg -simLearnType Full -verboseOutput -maxThreads 1 -bfile "'$snpFile'" -pheno "'$permPheno'" -mpheno "$pheno" -out "$batch"_perm_"'$chunk'" -numjobs 20 -thisjob "'$chunk'"
	done
	" >> $target

cd $oFolder ; $runSh $batch".fastlmm.sh" ; cd ..
#cd $oFolder ; sh $batch".fastlmm.sh" ; cd ..
done


iPhenotype=$bDir"/NewPhenotypeFile"
phenotype=$bDir/phenotype_res

head -1 $iPhenotype  > $bDir/tmp
nbCols=$(wc $bDir/tmp | awk '{print $2}')
nbGroups=$(expr $nbCols - 2)


for pheno in $(seq 1 $nbGroups)
do

tt=$(expr $pheno + 2)
batch=$(awk "{print \$$tt}" $bDir/tmp)
target=$oFolder"/"$batch".res.fastlmm.sh"; if [ -e $target ]; then rm $target; fi


echo "
	release=February2015

	bDir=/scratch2/vyp-scratch2/cian/UCLex_${release}/

	fastlmm='/share/apps/bahler/FaSTLMM.207/Bin/Linux_MKL/fastlmmc'
	snpFile=$bDir"allChr_snpStats_out"
	kinshipFile=$bDir"TechnicalKinship_Fastlmm"
	extract=$bDir"Clean_variants_Func"
	resPheno=$bDir"phenotype_res"
	for chunk in {1..20}
	do
	"'$fastlmm'"  -linreg -simLearnType Full -verboseOutput -maxThreads 1 -bfile "'$snpFile'" -pheno "'$resPheno'" -mpheno "$pheno" -out "$batch"_res_"'$chunk'" -numjobs 20 -thisjob "'$chunk'"
	done
	" >> $target

cd $oFolder ; $runSh $batch".res.fastlmm.sh" ; cd ..
#cd $oFolder ; sh $batch".fastlmm.sh" ; cd ..
done

