#!/bin/bash

ldak=/cluster/project8/vyp/cian/support/ldak/ldak

rootODir=${1-$rootODir}
release=${2-$release}
plink=/home/sejjcmu/bin/plink/plink
bDir=${rootODir}/UCLex_${release}/

## if [ ! -e ${bDir}/read_depth ]; then mkdir ${bDir}/read_depth; fi ## Not needed as I dont use read depth info right now and bim and fam are now made in first.step.R

## genotype to missing non missing
GenotypeMatrix=${rootODir}allChr_snpStats
missingNonMissing=${rootODir}Matrix.calls.Missing.NonMissing
ReadDepth=${rootODir}/read_depth/Depth_Matrix
phenFile=${rootODir}Phenotypes 
echo "Working with genotype matrix $GenotypeMatrix"

sed 's/0/2/g' $GenotypeMatrix".sp" | sed 's/1/2/g' | sed 's/NA/0/g' > ${missingNonMissing}.sp
ln -s ${rootODir}UCLex_${release}.bim ${ReadDepth}.bim
ln -s ${rootODir}UCLex_${release}.fam ${ReadDepth}.fam
ln -s ${rootODir}UCLex_${release}.bim ${missingNonMissing}.bim
ln -s ${rootODir}UCLex_${release}.fam ${missingNonMissing}.fam
$ldak --make-bed $missingNonMissing --sp $missingNonMissing


