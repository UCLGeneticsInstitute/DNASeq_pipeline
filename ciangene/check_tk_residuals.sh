#!/bin/bash

ldak=/cluster/project8/vyp/cian/support/ldak/ldak
#rootODir=/scratch2/vyp-scratch2/ciangene
rootODir=/scratch2/vyp-scratch2/cian/
release=February2015
rootODir=${1-$rootODir}
release=${2-$release}
R=/share/apps/R-3.1.0/bin/R
bDir=${rootODir}/UCLex_${release}/
genes=/SAN/biomed/biomed14/vyp-scratch/cian/LDAK/genesldak_ref.txt
kinship=$bDir"TechKin"
pKinship=$bDir"PopKin"
dKinship=$bDir"read_depth/Depthkin"
data=$bDir"allChr_snpStats_out"
phenotypes=$bDir"Phenotypes"
groups=$bDir"cohort.summary"

Names=$bDir"GroupNames"
awk '{print $4}' $groups > tmp
tail -n +2 "tmp" > $Names
rm tmp
nbGroups=$(wc -l $Names | awk {'print $1}') 

oDir=$bDir"KinshipDecomposition/"
if [ ! -e $oDir ]; then mkdir $oDir; fi

startPheno=1

for pheno in $(seq $startPheno $nbGroups)
do

	batch=$(sed -n $pheno'p' $Names); echo $batch is nb $pheno
	if (($pheno==$startPheno))
	then
		$ldak --reml $oDir$batch"_tech" --grm $kinship  --pheno $phenotypes --mpheno $pheno --eigen-save $oDir/techEigen
		$ldak --reml $oDir$batch"_geno" --grm $pKinship  --pheno $phenotypes --mpheno $pheno --eigen-save $oDir/popEigen
		$ldak --reml $oDir$batch"_tech" --grm $dKinship  --pheno $phenotypes --mpheno $pheno --eigen-save $oDir/depthEigen
	fi
		$ldak --reml $oDir$batch"_tech" --grm $kinship  --pheno $phenotypes --mpheno $pheno --eigen $oDir/techEigen	
		$ldak --reml $oDir$batch"_geno" --grm $pKinship  --pheno $phenotypes --mpheno $pheno --eigen $oDir/popEigen
		$ldak --reml $oDir$batch"_geno" --grm $dKinship  --pheno $phenotypes --mpheno $pheno --eigen $oDir/depthEigen

	if (($pheno==$startPheno))
	then
		awk '{ print $1, $1}' $oDir$batch"_tech.indi.res" > $bDir".NewPhenotypeFiletmp"
	fi
	
	cut -f5 $oDir$batch'_tech.indi.res' > .tmp
	tt=$(expr $pheno + 2)
	cat .tmp | cut -d ' ' -f $tt |  sed "1s/.*/$batch/"  > .tmp2
	paste  $bDir".NewPhenotypeFiletmp" .tmp2 > .tmp3
	mv .tmp3 $bDir".NewPhenotypeFiletmp" 
	

done

mv $bDir".NewPhenotypeFiletmp"  $bDir"NewPhenotypeFile"
tail -n +2 $bDir"NewPhenotypeFile" > $bDir"phenotype_res"

oFile=$oDir/plot.residuals.R
echo "dir<-'"$oDir"'" > $oFile
echo '
	files <- list.files(dir, pattern = "indi\\.res", full.names=T)

	pdf(paste0(dir, "/Residuals.pdf") )
		par(mfrow=c(2,2))  
		lapply(files, function(x)
		{
			name <- gsub(basename(x) , pattern = "\\.indi\\.res", replacement = "") 
			message(name)
			file <- read.table(x, header=T ) 
			# with(file, plot(Phenotype, Residual), xlab = name)
			plot(file$Phenotype, file$Residual, xlab = name) 
		}
		)
	dev.off() 

	' >> $oFile
$R CMD BATCH --no-save --no-restore $oFile



