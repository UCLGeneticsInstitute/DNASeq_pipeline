#!/bin/bash

ldak=/cluster/project8/vyp/cian/support/ldak/ldak
R=/share/apps/R/bin/R
rootODir=${1-$rootODir}
release=${2-$release}
bDir=${rootODir}/

data=/SAN/vyplab/UCLex/mainset_${release}/cian/allChr_snpStats_out
extract=../PCA/SNPs_for_pca
plink=/share/apps/genomics/plink-1.09beta/plink

UCLex_bed=${bDir}UCLex${release}_pca
OneKG_sp=${bDir}onekg_calls_for_uclex_snps
OneKG_bed=${bDir}OneKG

$ldak --make-bed $UCLex_bed --bfile $data --extract $extract
$ldak --make-bed $OneKG_bed --sp $OneKG_sp
$plink --noweb --bfile ${UCLex_bed}_out --bmerge ${OneKG_bed}_out.bed ${OneKG_bed}_out.bim ${OneKG_bed}_out.fam --make-bed --out ${bDir}UCLex_OneKG_merged --remove ${bDir}onekg.samples.in.ucl
$ldak --calc-kins-direct ${bDir}UCLex_OneKG_merged_kin --bfile ${bDir}UCLex_OneKG_merged --ignore-weights YES
$ldak --pca ${bDir}UCLex${release}_OneKG_merged_pca --bfile ${bDir}UCLex_OneKG_merged --grm ${bDir}UCLex_OneKG_merged_kin



echo '5-getPCAsnps' >> $bDir/Check
