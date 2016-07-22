#!/bin/bash

ldak=/cluster/project8/vyp/cian/support/ldak/ldak.4.8
R=/share/apps/R-3.1.0/bin/R
#rootODir=/scratch2/vyp-scratch2/ciangene
rootODir=/scratch2/vyp-scratch2/cian/
release=June2015
rootODir=${1-$rootODir}
release=${2-$release}
bDir=${rootODir}/UCLex_${release}/


data=/scratch2/vyp-scratch2/cian/UCLex_${release}/allChr_snpStats_out
ldak=/cluster/project8/vyp/cian/support/ldak/ldak
extract=/cluster/project8/vyp/cian/data/UCLex/ciangene/scripts/PCA/SNPs_for_pca
plink=/share/apps/genomics/plink-1.07-x86_64/plink

UCLex_bed=${bDir}UCLex${release}_pca
OneKG_sp=${bDir}onekg_calls_for_uclex_snps
OneKG_bed=${bDir}OneKG

$plink --noweb --bfile ${UCLex_bed}_out --bmerge ${OneKG_bed}_out.bed ${OneKG_bed}_out.bim ${OneKG_bed}_out.fam --make-bed --out ${bDir}UCLex_OneKG_merged  --remove ${bDir}onekg.samples.in.ucl
$ldak --calc-kins-direct ${bDir}UCLex_OneKG_merged_kin --bfile ${bDir}UCLex_OneKG_merged --ignore-weights YES
$ldak --pca ${bDir}UCLex${release}_OneKG_merged_pca --bfile ${bDir}UCLex_OneKG_merged --grm ${bDir}UCLex_OneKG_merged_kin
