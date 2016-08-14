#/bin/sh

rootODir=/SAN/vyplab/UCLex/mainset_July2016/cian/
release=July2016

rootODir=${1-$rootODir}
release=${2-$release}

runSh='sh /cluster/project8/vyp/cian/scripts/bash/runBashCluster.sh'
Groups=${rootODir}GroupNames
nbGroups=$(wc -l  $Groups | awk '{print $1}') 
templateScript="/SAN/vyplab/UCLex/scripts/DNASeq_pipeline/ciangene/Gene_Based_Tests/LDAK/ldak_base.sh"
variants=${rootODir}filtered_snps_skat
phenotype=${rootODir}Phenotypes

oFolder=${rootODir}LDAK
if [ ! -e $oFolder ]; then mkdir $oFolder; fi
cd $oFolder 

## Chosen Kinship 
##  Model 					MeanVarianceExplained   MAF 	CallRate	FILTER
##  maf_0.001_callRate_0.5_clean.progress       0.7577674 		1e-03 	0.50000		TRUE

nbGroups=106
for pheno in $(seq 73 $nbGroups) 
do
	batch=$(sed -n $pheno'p' $Groups)
	for maf in {0.00001,0.001,,0.01}
	do

		for missing in {0.00001,0.5,0.9}
		do
				oFile=${batch}_${maf}_${missing}_base.LDAK.sh
				target=$oFolder/$oFile
				echo "minMaf='$maf' ; "'minMaf=$(echo $minMaf/10|bc -l)'"; "'minMaf=$(echo $minMaf | sed 's/0*$//')'" " > $target	
				echo "missing='$missing' ; "'missing=$(echo $missing/10|bc -l)'"; "'missing=$(echo $missing | sed 's/0*$//')'"" >> $target	
				echo "phenotypes='$phenotype' ; mPheno=$pheno ; pheno='$batch'"  >> $target
				echo "extract='$variants' ; role='Func'" >> $target
				cat $templateScript >> $target
				$runSh $oFile
		done
	done
done

