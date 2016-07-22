#/bin/sh

#rootODir=/scratch2/vyp-scratch2/ciangene
rootODir=/scratch2/vyp-scratch2/cian/
release=February2015
rootODir=${1-$rootODir}
release=${2-$release}
bDir=${rootODir}/UCLex_${release}/

runSh='sh /cluster/project8/vyp/cian/scripts/bash/runBashCluster.sh'
Groups=$bDir"GroupNames"
nbGroups=$(wc -l  $Groups | awk '{print $1}') 
#templateScript="/cluster/project8/vyp/cian/data/UCLex/UCLex_August/Scripts/ciangene/scripts/LDAK/DepthKin_NormalPheno_template.sh"
templateScript=/cluster/project8/vyp/cian/data/UCLex/ciangene/scripts/LDAK/DepthKin_NormalPheno_template.sh 
variants=$bDir"/Clean_variants_func_rare"
phenotype=$bDir"/Phenotypes"

oFolder=$bDir"/LDAK_gene_tests_all_phenos/"
if [ ! -e $oFolder ]; then mkdir $oFolder; fi

for pheno in $(seq 80 $nbGroups) 
#for pheno in {79,90,91,102} # skipping the first few random phenos
do
	batch=$(sed -n $pheno'p' $Groups)
	for maf in {0.0000001,0.001}
	do

		for missing in {0.000001,9}
		do
			oFile=${batch}_${maf}_${missing}.LDAK.sh
			target=$oFolder/$oFile
			echo "minMaf='$maf' ; "'minMaf=$(echo $minMaf/10|bc -l)'"; "'minMaf=$(echo $minMaf | sed 's/0*$//')'" " > $target	
			echo "missing='$missing' ; "'missing=$(echo $missing/10|bc -l)'"; "'missing=$(echo $missing | sed 's/0*$//')'"" >> $target	
			echo "phenotypes='$phenotype' ; mPheno=$pheno ; pheno='$batch'"  >> $target
			echo "extract='$variants' ; role='Func'" >> $target
			cat $templateScript >> $target
			cd $oFolder ; $runSh $oFile ; cd ..
			#sh $target
			done
	done
done

