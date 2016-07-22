#!/bin/bash
shopt -s expand_aliases
source ~/.bashrc

release=October2014
snps="/cluster/project8/vyp/cian/data/UCLex/UCLex_${release}/All_phenotypes/LDAK/LDAK_gene_tests_all_phenos/no_kin_maf_0.01_Func_ALevine/gene_preds.txt"
genoFile="/scratch2/vyp-scratch2/cian/UCLex_${release}/Genotype_Matrix.sp_out"
techFile="/scratch2/vyp-scratch2/cian/UCLex_${release}/missingNonmissingPlink_out"

genoOutFile="/scratch2/vyp-scratch2/cian/UCLex_${release}/Genotype_Matrix_filt"
techOutFile="/scratch2/vyp-scratch2/cian/UCLex_${release}/missingNonmissing_filt"
plink --noweb --bfile $genoFile --extract $snps --make-bed --out $genoOutFile
plink --noweb --bfile $techFile --extract $snps --make-bed --out $techOutFile

techKinOut="/scratch2/vyp-scratch2/cian/UCLex_${release}/TechKin_filt"
techPCsOut="/scratch2/vyp-scratch2/cian/UCLex_${release}/TechPCs_filt"
ldak --calc-kins-direct $techKinOut --bfile $techOutFile --ignore-weights YES --kinship-raw YES
ldak --pca $techPCsOut --grm $techKinOut --bfile  $techOutFile --ignore-weights YES
