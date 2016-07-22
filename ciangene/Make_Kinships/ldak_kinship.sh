ldak=/cluster/project8/vyp/cian/support/ldak/ldak
## Test variance explained by multiple kinships
$ldak --reml test --mgrm list --pheno /cluster/project8/vyp/cian/data/UCLex/UCLex_February2015/LambiaseTechKin/lambiase.pheno
## Find snps that drive kinship signal 
$ldak --calc-pca-loads  out --grm TechKin --bfile Matrix.calls.Missing.NonMissing_out --pcastem  TechKin
