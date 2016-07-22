#!/bin/bash

ldak=/cluster/project8/vyp/cian/support/ldak/ldak
#rootODir=/scratch2/vyp-scratch2/ciangene
rootODir=/scratch2/vyp-scratch2/cian/
release=June2015
rootODir=${1-$rootODir}
release=${2-$release}
plink=/share/apps/genomics/plink-1.07-x86_64/plink
bDir=${rootODir}/UCLex_${release}/

## if [ ! -e ${bDir}/read_depth ]; then mkdir ${bDir}/read_depth; fi ## Not needed as I dont use read depth info right now and bim and fam are now made in first.step.R

## genotype to missing non missing
GenotypeMatrix=${bDir}allChr_snpStats
missingNonMissing=${bDir}Matrix.calls.Missing.NonMissing
ReadDepth=${bDir}/read_depth/Depth_Matrix
phenFile=${bDir}Phenotypes 
echo "Working with genotype matrix $GenotypeMatrix"

sed 's/0/2/g' $GenotypeMatrix".sp" | sed 's/1/2/g' | sed 's/NA/0/g' > ${missingNonMissing}.sp
ln -s ${bDir}UCLex_${release}.bim ${ReadDepth}.bim
ln -s ${bDir}UCLex_${release}.fam ${ReadDepth}.fam
ln -s ${bDir}UCLex_${release}.bim ${missingNonMissing}.bim
ln -s ${bDir}UCLex_${release}.fam ${missingNonMissing}.fam
$ldak --make-bed $missingNonMissing --sp $missingNonMissing


