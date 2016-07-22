ldak=/cluster/project8/vyp/cian/support/ldak/ldak

rootODir=/scratch2/vyp-scratch2/ciangene
release=February2015


## got to update this
	# DO i still want to bother calculating genotype weights if im not using them 
	# no need for read depth kin
	## use calc-kins-direct instead of cutting and joining. probably slower but more convenient. 



rootODir=$1
release=$2


bDir=${rootODir}/UCLex_${release}

missingNonMissing=$bDir"/Matrix.calls.Missing.NonMissing.sp"



###################################################################################################################################################
##  I should calculate 3 sets of SNP weights:
	# 1. Genotype based weights for the actual gene based tests
	# 2. MissingNonMissing weights to be used when making the technical kinship
	# 3. ReadDepth weights to be used when making the DepthKin


## Some basic parameters: 
minObs=0.9 			## SNP needs to be present in 90% samples to be included. 
minMaf=0.000001			## SNP with MAF >= this are retained
minVar=0.001			## SNP with variance >= this are retained? 
maxTime=500			## Nb minutes calculation allowed run for. 



################################################################################################################################################
### Lets get started, calculate genotype weights
################################################################################################################################################
data=$bDir"Genotype_Matrix_out" 
weights=$bDir"Genotype_weights"

rm  $data".bim" ; ln -s $bDir"read_depth/Depth_Matrix.bim" $data".bim"
rm  $data".fam" ; ln -s $bDir"read_depth/Depth_Matrix.fam" $data".fam"

## cut data into sections and calculate SNP weights for each. 
$ldak --cut-weights  $weights  --bfile $data
regions=$(tail -1 $weights/section_details.txt | awk '{print $1}')
mkdir $bDir"Scripts"
for region in $(seq 1 $regions)
	do
	oFile=$bDir"Scripts/calc_weights_"$region".sh"
	echo "$ldak --calc-weights $weights --bfile $data --section $region --maxtime $maxTime --minmaf $minMaf --minvar $minVar --minobs $minObs " > $oFile
	runSh $oFile
done



## Merge all the weights together
oFile="./scripts/merge_weights.sh"
echo "
$ldak --join-weights $weights --bfile $data 
" > $oFile
runSh $oFile


## re calculate the weights, recommended by Doug
$ldak --cut-weights $weights --bfile $data  --weights $weights/weightsALL 

regions=$(tail -1 $weights/re-section_details.txt | awk '{print $1}')
for region in $(seq 1 $regions)
	do
	oFile=$bDir"Scripts/re_calc_weights_"$region".sh"
	echo "$ldak --calc-weights $weights --bfile $data --section $region --maxtime $maxTime --minmaf $minMaf --minvar $minVar --minobs $minObs --weights $weights/weightsALL " > $oFile
	runSh $oFile
done
$ldak --join-weights $weights --bfile $data --weights $weights/weightsALL 




####################################################################################################################
#####   Calculate tech snp weights and kinship
###################################################################################################################

data=$bDir"Matrix.calls.Missing.NonMissing" 
weights=$bDir"TechKin_weights"

rm  $data".bim" ; ln -s $bDir"read_depth/Depth_Matrix.bim" $data".bim"
rm  $data".fam" ; ln -s $bDir"read_depth/Depth_Matrix.fam" $data".fam"


## cut data into sections and calculate SNP weights for each. 
$ldak --cut-weights  $weights  --sp $data
regions=$(tail -1 $weights/section_details.txt | awk '{print $1}')
for region in $(seq 1 $regions)
	do
	oFile=$bDir"Scripts/calc_weights_"$region".sh"
	echo "$ldak --calc-weights $weights --sp $data --section $region --maxtime $maxTime  --minmaf $minMaf --minvar $minVar --minobs $minObs " > $oFile
	runSh $oFile
done


## Merge all the weights together
oFile="./scripts/merge_weights.sh"
echo "
$ldak --join-weights $weights --sp $data 
" > $oFile
runSh $oFile



$ldak --cut-weights $weights --sp $data  --weights $weights/weightsALL 

## We've cut into sections, now re-calculating the weights for each section
regions=$(tail -1 $weights/re-section_details.txt | awk '{print $1}')
for region in $(seq 1 $regions)
	do
	oFile=$bDir"Scripts/re_calc_weights_"$region".sh"
	echo "$ldak --calc-weights $weights --sp $data --section $region --maxtime $maxTime --minmaf $minMaf --minvar $minVar --minobs $minObs --weights $weights/weightsALL " > $oFile
	runSh $oFile
 done


$ldak --join-weights $weights --sp $data --weights $weights/weightsALL 


## Make Technical Kinship  matrix
$ldak --cut-kins $weights --sp $data

regions=$(tail -1 $weights/partition_details.txt | awk '{print $1}') 

for region in $(seq 1 $regions)
	do
	oFile=$bDir"Scripts/kin_calc_"$region".sh"
	echo "$ldak --calc-kins $weights  --partition $region --sp $data --weights $weights/re-weightsALL --power 0 " > $oFile
	runSh $oFile
done

oFile="$oFolder/join_regions.sh"
echo "$ldak --join-kins $weights --kinship-matrix YES" > $oFile
runSh $oFile



####################################################################################################################
#####   Calculate Read Depth SNP weights
###################################################################################################################

data=$bDir"read_depth/Depth_Matrix"
weights=$bDir"Depth_weights"

## cut data into sections and calculate SNP weights for each. 
$ldak --cut-weights  $weights  --sp $data
regions=$(tail -1 $weights/section_details.txt | awk '{print $1}')
mkdir $bDir"Scripts"
for region in $(seq 1 $regions)
	do
	oFile=$bDir"Scripts/calc_weights_"$region".sh"
	echo "$ldak --calc-weights $weights --sp $data --section $region --minmaf $minMaf --maxtime $maxTime --minvar $minVar --minobs $minObs " > $oFile
	runSh $oFile
done

oFile="./scripts/merge_weights.sh"
echo "
$ldak --join-weights $weights --sp $data 
" > $oFile
runSh $oFile


$ldak --cut-weights $weights --sp $data  --weights $weights/weightsALL 

## We've cut into sections, now re-calculating the weights for each section
regions=$(tail -1 $weights/re-section_details.txt | awk '{print $1}')
for region in $(seq 1 $regions)
	do
	oFile=$bDir"Scripts/re_calc_weights_"$region".sh"
	echo "$ldak --calc-weights $weights --sp $data --section $region --maxtime $maxTime --minmaf $minMaf --minvar $minVar --minobs $minObs --weights $weights/weightsALL " > $oFile
	runSh $oFile
done
$ldak --join-weights $weights --sp $data --weights $weights/weightsALL 



############# Calculate Read Depth Kinship matrix
$ldak --cut-kins $weights --sp $data

regions=$(tail -1 $weights/partition_details.txt | awk '{print $1}') 

for region in $(seq 1 $regions)
	do
	oFile=$bDir"Scripts/kin_calc_"$region".sh"
	echo "$ldak --calc-kins $weights  --partition $region --sp $data --weights $weights/re-weightsALL --power 0 " > $oFile
	runSh $oFile
done

oFile="$oFolder/join_regions.sh"
echo "
$ldak --join-kins $weights --kinship-matrix YES
" > $oFile
runSh $oFile






