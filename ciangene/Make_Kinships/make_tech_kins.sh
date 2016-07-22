#/bin/sh
shopt -s expand_aliases
source ~/.bashrc

# The genotype data has already been converted to the missingNonMissing format. Here, I extract the samples and SNPs i want and then calculate both the Technical Kinship matrices and the Principal Components. 

release=October2014
datDir="/scratch2/vyp-scratch2/cian/UCLex_${release}"
data=$datDir"/missingNonmissingPlink_out"
bDir="/cluster/project8/vyp/cian/data/UCLex/UCLex_${release}/Gene_tests/Clean_Up"
AllKeep='/cluster/project8/vyp/cian/data/UCLex/UCLex_${release}/Lambiase_case_control/whole_exome/support/Lambiase_joint_ivf_sad_keep'
oDir='/scratch2/vyp-scratch2/cian/UCLex_${release}/Lambiase_case_control'



# ivf targeted snps
plink --noweb --bfile $data --keep $bDir/ivfKeep --make-bed --extract SNPS.in.target.genes --range --out $oDir/IVFvsUCLexTargetedGenes_missingNonMissing 
ldak --calc-kins-direct $oDir/IVFvsUCLexTargetedGenes_missingNonMissing_tech_kin  --bfile $oDir/IVFvsUCLexTargetedGenes_missingNonMissing --ignore-weights YES --kinship-raw YES
ldak --pca $oDir/IVFvsUCLexTargetedGenes_missingNonMissing_tech_pcs --bfile $oDir/IVFvsUCLexTargetedGenes_missingNonMissing --ignore-weights YES --grm $oDir/IVFvsUCLexTargetedGenes_missingNonMissing_tech_kin 

# lambiase joint all snps
plink --noweb --bfile $data --keep $AllKeep --make-bed --out $oDir/LambiaseVsUCLex_allSnps_missingNonMissing
ldak --calc-kins-direct $oDir/LambiaseVsUCLex_allSnps_missingNonMissing_tech_kin --bfile $oDir/LambiaseVsUCLex_allSnps_missingNonMissing --ignore-weights YES --kinship-raw YES
ldak --pca $datDir/IVFvsUCLex_tech_pcs --bfile $datDir/IVFvsUCLex --ignore-weights YES --grm $datDir/IVFvsUCLex_tech_kin
