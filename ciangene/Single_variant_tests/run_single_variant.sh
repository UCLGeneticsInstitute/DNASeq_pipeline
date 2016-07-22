#!/bin/bash
shopt -s expand_aliases
source ~/.bashrc

release=February2015
bDir=/scratch2/vyp-scratch2/cian/UCLex_${release}/
Pheno=$bDir"Phenotypes"
Groups=$bDir"GroupNames"
phenos=$(wc -l  $Groups | awk '{print $1}') 

mkdir bin

for fil in {1..23}
do
	for phen in  $(seq 1 $phenos)
	do
		oFile=$(sed -n $phen'p' Groups)
		ivfOut=bin/UCLex_fisher_${oFile}_${fil}.R
		echo "fil <- $fil"  > $ivfOut
		echo "pheno  <- '$Pheno'" >> $ivfOut
		echo "phenoType  <- $phen" >> $ivfOut
		echo "oFile  <- '$oFile'" >> $ivfOut
		cat single_variant_fisher.R >> $ivfOut
		target=$(basename $ivfOut)
		cd bin; runR $target; cd ..
	done
done
