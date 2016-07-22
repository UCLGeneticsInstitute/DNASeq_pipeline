ciangene
========

[![Join the chat at https://gitter.im/CianMurphy/ciangene](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/CianMurphy/ciangene?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
# Order of scripts: 

* first.step.R

		Takes the snpStats files from the data release and turn into a numeric matrix
		Establishes UCLex External Controls (Uex) and calculates their summary stats
		Creates ReadDepth.sp

* variant_filtering/qc.sh 
		
		filters by Uex missingness and MAF
		Makes geno bed file

* variant_filtering/filter_snps.R

		filters by function - gene tests

* PCA 

		Identify UCLex sample ancestry by comparing to 1000G. Remove non Caucasians from further analyses. 

* new_make_phenotype_file.R

		Creates a large matrix of phenotypes, grouped by names basically.
		Fix sample groupings
		For each phenotype, create external control set. Use for QC and for MAF filter. 
		For each phenotype, remove confounding controls. 

* Case control missingness

		For every selected phenotype, identify variants that have a significantly different missingness rate between cases and controls. Useful for subsequent QC. 

* convert_genotype_to_missingNonMissing.sh

		turns NA/non NA into 0-1 for LDAK/plink use
		Makes TK bed 

* make_kinships_new.sh

		Makes technical kinship (TK), Read Depth (RD) and geno kinship matrix
		calculate and plot PCs

* check_tk_residuals.sh
		
		Checks how much kinships explain for all phenos. 
		Creates new pheno files, using residuals from (B) as phenos

* LDAK/run_ldak_on_all_phenos.sh & run_ldak_on_all_phenos_res.sh

		run_ldak_on_all_phenos.sh - gene tests - all UCLex phenos, with and without kinship
		run_ldak_on_all_phenos_res.sh - gene tests using new pheno file with TK residuals. 		

* plot_gene_test.R
		
		plot gene based results from LDAK. 
	
* plink_single_variant_tests.sh

		single variant tests - Fisher, logistic regression with and without techPC covariates 

* run_fastLMM_on_all_phenotypes.sh 
	
		first, run prepare_kinship_matrix_for_fastLMM.R

		single variant tests - FastLMM - Linear Regression / Mixed model with and without techKin


* merge_gene_tests.sh
	
		Doesn't just merge. Also uses permuted pvalue in regress1..regressN to define null distribution (permuted pheno with kin) to calc final pvalue in regressALL. Script hangs when run on just one chunk/gene (eg as per power studies). 

* plot_singleVariant_results.R

		Merges single variant tests and plots some of them. 

* annotate_qqplot.R 

		Makes qqPlots of snpStats style but allows you to label points too. 
