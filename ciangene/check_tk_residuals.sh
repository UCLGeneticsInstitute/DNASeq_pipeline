#!/bin/bash

ldak=/cluster/project8/vyp/cian/support/ldak/ldak


rootODir=${1-$rootODir}
release=${2-$release}
R=/share/apps/R/bin/R

genes=/SAN/biomed/biomed14/vyp-scratch/cian/LDAK/genesldak_ref.txt
kinship=${rootODir}TechKin
pKinship=${rootODir}PopKin
dKinship=${rootODir}DepthKin
data=${rootODir}allChr_snpStats_out
phenotypes=${rootODir}Phenotypes
groups=${rootODir}cohort.summary

Names=${rootODir}GroupNames
awk '{print $4}' $groups > tmp
tail -n +2 "tmp" > $Names
rm tmp
nbGroups=$(wc -l $Names | awk {'print $1}') 

oDir=${rootODir}KinshipDecomposition/
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
		awk '{ print $1, $1}' $oDir$batch"_tech.indi.res" > ${rootODir}.NewPhenotypeFiletmp
	fi
	
	cut -f5 $oDir$batch'_tech.indi.res' > .tmp
	tt=$(expr $pheno + 2)
	cat .tmp | cut -d ' ' -f $tt |  sed "1s/.*/$batch/"  > .tmp2
	paste  ${rootODir}.NewPhenotypeFiletmp .tmp2 > .tmp3
	mv .tmp3 ${rootODir}.NewPhenotypeFiletmp 
	

done

mv ${rootODir}.NewPhenotypeFiletmp ${rootODir}NewPhenotypeFile
tail -n +2 ${rootODir}NewPhenotypeFile > ${rootODir}phenotype_res

oFile=$oDir/plot.residuals.R
pDir=${rootODir}plots/
echo "pdir<-'"$pDir"'" > $oFile
echo "dir<-'"$oDir"'" >> $oFile
echo '
	files <- list.files(dir, pattern = "indi\\.res", full.names=T)

	pdf(paste0(pdir, "/Residuals.pdf") )
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
