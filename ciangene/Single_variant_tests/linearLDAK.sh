#rootODir=/scratch2/vyp-scratch2/ciangene
rootODir=/scratch2/vyp-scratch2/cian/
release=July2015
rootODir=${1-$rootODir}
release=${2-$release}
bDir=${rootODir}/UCLex_${release}/LDAK_gene_tests_all_phenos/
runSh='sh /cluster/project8/vyp/cian/scripts/bash/runBashCluster.sh'
phenoFile=/scratch2/vyp-scratch2/cian/UCLex_${release}/Clean_pheno_subset_tk_depth_residuals
CleanPheno=/scratch2/vyp-scratch2/cian/UCLex_${release}/Clean_pheno_subset

#### Data
ldak=/cluster/project8/vyp/cian/support/ldak/ldak
data=$rootODir/UCLex_${release}/allChr_snpStats_out
kinship=/cluster/project8/vyp/cian/data/UCLex/ciangene/scripts/TK_maf_0.1_callRate_0.4_clean
depthKin=/cluster/project8/vyp/cian/data/UCLex/ciangene/scripts/Depth_maf_0.1_callRate_0.00001_clean
genes=/SAN/biomed/biomed14/vyp-scratch/cian/LDAK/genesldak_ref.txt

####### Parameters
#minMaf=0.000001
#missing=0.7
minVar=1e-4
maxmaf=0.5
power=0 # rare variants arent upweighed with power = 0
minWeight=3 # genes need min this many variants to be included
overlap=NO # variants counted once per region, no overlapping transcripts
ignoreWeights=YES ## if YES, variants arent weighted for LD correction
partition=1
#######
kinships=/cluster/project8/vyp/cian/data/UCLex/ciangene/scripts/KinshipList

outDir=${rootODir}/UCLex_${release}/LDAK_Single_Variant_Tests/
if [ ! -e $outDir ]; then mkdir $outDir; fi

for pheno in {1..16}
do
	oFile=${outDir}test${pheno}_base.sh
	echo "
	$ldak --linear ${outDir}pheno${pheno}_base --pheno $CleanPheno --bfile $data --proximal NO --prox-buffer 0 --mpheno $pheno 
	" > $oFile
	$runSh $oFile

	oFile=${outDir}test${pheno}_base_kin.sh
	echo "
	$ldak --linear ${outDir}pheno${pheno}_tk --pheno  $CleanPheno --bfile $data --grm $kinship --proximal NO --prox-buffer 0 --mpheno $pheno 
	" > $oFile
	$runSh $oFile

	oFile=${outDir}test${pheno}_residuals.sh
	echo "
	$ldak --linear ${outDir}pheno${pheno}_residuals --pheno  $phenoFile --bfile $data --proximal NO --prox-buffer 0 --mpheno $pheno 
	" > $oFile
	$runSh $oFile

	oFile=${outDir}test${pheno}_base_dual_kin.sh
	echo "
	$ldak --linear ${outDir}pheno${pheno}_dual_kin --pheno  $CleanPheno --bfile $data --mgrm $kinships --proximal NO --prox-buffer 0 --mpheno $pheno 
	" > $oFile
	$runSh $oFile

	oFile=${outDir}test${pheno}_perm.sh
	echo "
	$ldak --linear ${outDir}pheno${pheno}_perm --pheno  $CleanPheno --bfile $data --proximal NO --prox-buffer 0 --mpheno $pheno --permute YES
	" > $oFile
	$runSh $oFile

	oFile=${outDir}test${pheno}_rd.sh
	echo "
	$ldak --linear ${outDir}pheno${pheno}_rd --pheno  $CleanPheno --bfile $data --proximal NO --prox-buffer 0 --mpheno $pheno --grm $depthKin
	" > $oFile
	$runSh $oFile
done

