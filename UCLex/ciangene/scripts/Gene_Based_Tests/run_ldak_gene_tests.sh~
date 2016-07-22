#!/bin/bash
shopt -s expand_aliases
source ~/.bashrc


Out=IVF_lof_LDAK.sh
echo "phenotypes='/cluster/project8/vyp/cian/data/UCLex/UCLex_October2014/Lambiase_case_control/support/IVF.pheno' ; pheno='IVFvsUCL'"  > $Out
echo "kinship='/scratch2/vyp-scratch2/cian/UCLex_October2014/Lambiase_case_control/IVF_all_kin'; kin='Tech'" >> $Out
echo "extract='/cluster/project8/vyp/cian/data/UCLex/UCLex_October2014/Lambiase_case_control/whole_exome/gene_based/SNPs.lof_extract' ; role='Lof'" >> $Out
cat ldak_gene_tests_template.sh >> $Out
sh $Out

Out=IVF_func_LDAK.sh
echo "phenotypes='/cluster/project8/vyp/cian/data/UCLex/UCLex_October2014/Lambiase_case_control/support/IVF.pheno' ; pheno='IVFvsUCL'"  > $Out
echo "kinship='/scratch2/vyp-scratch2/cian/UCLex_October2014/Lambiase_case_control/IVF_all_kin'; kin='Tech'" >> $Out
echo "extract='/cluster/project8/vyp/cian/data/UCLex/UCLex_October2014/Lambiase_case_control/whole_exome/gene_based/SNPs.func_extract' ; role='Func'" >> $Out
cat ldak_gene_tests_template.sh >> $Out
sh $Out


Out=Lambiase_lof_LDAK.sh
echo "phenotypes='/cluster/project8/vyp/cian/data/UCLex/UCLex_October2014/Lambiase_case_control/support/LambiaseVsUCLex.ex.ctrls.removed.pheno' ; pheno='LambiasevsUCL'"  > $Out
echo "kinship='/scratch2/vyp-scratch2/cian/UCLex_October2014/Lambiase_case_control/LambiaseVsUCLex_allSnps_missingNonMissing_tech_kin'; kin='Tech'" >> $Out
echo "extract='/cluster/project8/vyp/cian/data/UCLex/UCLex_October2014/Lambiase_case_control/whole_exome/gene_based/SNPs.lof_extract' ; role='lof'" >> $Out
cat ldak_gene_tests_template.sh >> $Out
sh $Out

Out=Lambiase_func_LDAK.sh
echo "phenotypes='/cluster/project8/vyp/cian/data/UCLex/UCLex_October2014/Lambiase_case_control/support/LambiaseVsUCLex.ex.ctrls.removed.pheno' ; pheno='LambiasevsUCL'"  > $Out
echo "kinship='/scratch2/vyp-scratch2/cian/UCLex_October2014/Lambiase_case_control/LambiaseVsUCLex_allSnps_missingNonMissing_tech_kin'; kin='Tech'" >> $Out
echo "extract='/cluster/project8/vyp/cian/data/UCLex/UCLex_October2014/Lambiase_case_control/whole_exome/gene_based/SNPs.func_extract' ; role='Func'" >> $Out
cat ldak_gene_tests_template.sh >> $Out
sh $Out

