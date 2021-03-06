#!/bin/bash

ldak=/home/sejjcmu/bin/ldak/ldak.4.9.fast
#ldak=/scratch2/vyp-scratch2/cian/ldak.5.98
R=/share/apps/R/bin/R
#rootODir=/scratch2/vyp-scratch2/ciangene


rootODir=${1-$rootODir}
release=${2-$release}
#bDir=${rootODir}/UCLex_${release}/

missingNonMissing=${rootODir}Matrix.calls.Missing.NonMissing_out
techOut=${rootODir}TechKin
Depth=${rootODir}/read_depth/Depth_Matrix
DepthOut=${rootODir}DepthKin
extract=${rootODir}Clean_variants_func_rare

## Some basic parameters: 
minObs=0.7	## SNP needs to be present in 90% samples to be included. 
minMaf=0.001			## SNP with MAF >= this are retained
maxMaf=0.5				## SNP with MAF <= this are retained # for techKIN maf is missingness rate. 
minVar=0.0000001			## SNP with variance >= this are retained? 
## maxTime=500			## Nb minutes calculation allowed run for. 


######### Tech Kin
$ldak --calc-kins-direct $techOut --bfile $missingNonMissing --ignore-weights YES --kinship-raw YES \
 --minmaf $minMaf --maxmaf $maxMaf --minvar $minVar --minobs $minObs --extract $extract 

$ldak --pca ${rootODir}TechPCs --grm $techOut 


oFile=${rootODir}plot.techpca.R
echo "dir<-'"${rootODir}"'" > $oFile
echo '
	file <- read.table(paste0(dir, "TechPCs.vect"), header=F) 
	groups <- read.table(paste0(dir, "Sample.cohort"), header=F)
	uniq.groups <- unique(groups[,2])
	nb.groups <- length(uniq.groups)

	buffer <- 0.01
	xmin <- min(file[,3]) - buffer
	xmax <- max(file[,3]) + buffer
	ymin <- min(file[,4]) - buffer
	ymax <- max(file[,4]) + buffer

	png(paste0(dir, "/plots/TechPCA.png") ) 
	for(i in 1:nb.groups)
#	for(i in 72:nb.groups)
	{
		hit <- which(groups[,2] == uniq.groups[i])
		if(i==1)
		{
		plot(file[hit,3], file[hit,4], 
			xlab = "PC1", ylab = "PC2", 
			xlim=c(xmin, xmax), 
			ylim=c(ymin, ymax), 
			main = paste("TechPCA", date()) ,
			col=i
			) 
		} else
		{
			points(file[hit,3], file[hit,4], col=i)
		}
	}
	dev.off()
	png(paste0(dir, "/plots/TechPCA.var.png") ) 
	file <- read.table(paste0(dir,"TechPCs.values"), header=F)
	file$sd <- file[,1]^.5
	file$var <- file[,1]^2 / sum(file[,1]^2)
	plot(file$var*100, xlab = "Technical Principal Components" , ylab = "Proportion Variance Explained (%)", main = "Variance explained by each Principal Component") 
	dev.off() 

	vec<-read.table(paste0(dir, "TechPCs.vect"),header=T,sep=" ")
	miss<-read.table(paste0(dir, "gstats.imiss"),header=T,sep="\t",row.names=NULL)
	plot.data<-data.frame(PC1=vec[,3],PC2=vec[,4],CallRate=miss$F_MISS)
	png(paste0(dir, "TechPCS_and_CallRate.png")) 
	print(qplot(PC1,PC2,color=CallRate,data=plot.data,main="Technical PCs distinguish samples based on Missingness"))
	dev.off()

	' >> $oFile
$R CMD BATCH --no-save --no-restore $oFile



######### Pop Kin
missingNonMissing=${rootODir}allChr_snpStats_out
techOut=${rootODir}PopKin

$ldak --calc-kins-direct $techOut --bfile $missingNonMissing --ignore-weights YES --kinship-raw YES \
--minmaf $minMaf --maxmaf $maxMaf --minvar $minVar --minobs $minObs  --extract $extract 
$ldak --pca ${rootODir}popPCs --grm $techOut  


oFile=${rootODir}plot.popPca.R
echo "dir<-'"${rootODir}"'" > $oFile
echo '
	file <- read.table(paste0(dir, "popPCs.vect"), header=F) 
	groups <- read.table(paste0(dir, "Sample.cohort"), header=F)
	uniq.groups <- unique(groups[,2])
	nb.groups <- length(uniq.groups)

	buffer <- 0.01
	xmin <- min(file[,3]) - buffer
	xmax <- max(file[,3]) + buffer
	ymin <- min(file[,4]) - buffer
	ymax <- max(file[,4]) + buffer

	png(paste0(dir, "/plots/popPCA.png") ) 
	for(i in 1:nb.groups)
#	for(i in 72:nb.groups)
	{
		hit <- which(groups[,2] == uniq.groups[i])
		if(i==1)
		{
		plot(file[hit,3], file[hit,4], 
			xlab = "PC1", ylab = "PC2", 
			xlim=c(xmin, xmax), 
			ylim=c(ymin, ymax), 
			main = paste("popPCA", date()) ,
			col=i
			) 
		} else
		{
			points(file[hit,3], file[hit,4], col=i)
		}
	}
	dev.off()
	png(paste0(dir, "/plots/PopPCA.var.png") ) 
	file <- read.table(paste0(dir,"popPCs.values"), header=F)
	# file$sd <- file[,1]^.5
	file$var <- file[,1]^2 / sum(file[,1]^2)
	plot(file$var*100, xlab = "Genotype Principal Components" , ylab = "Proportion Variance Explained (%)", main = "Variance explained by each Principal Component") 
	dev.off() 

	' >> $oFile
$R CMD BATCH --no-save --no-restore $oFile


######### Depth Kin

## Some basic parameters: 
minObs=0 ## SNP needs to be present in 90% samples to be included. 
minMaf=0			## SNP with MAF >= this are retained
maxMaf=0.5				## SNP with MAF <= this are retained # for techKIN maf is missingness rate. 
minVar=0.0000001			## SNP with variance >= this are retained? 
## maxTime=500			## Nb minutes calculation allowed run for. 


DepthOut=${rootODir}DepthKin
$ldak --calc-kins-direct $DepthOut --sp $Depth --ignore-weights YES --kinship-raw YES \
--minmaf $minMaf --minvar $minVar --minobs $minObs  --extract $extract 
$ldak --pca ${rootODir}DepthPCs --grm $DepthOut 


oFile=${rootODir}plot.Depthpca.R
echo "dir<-'"${rootODir}"'" > $oFile
echo '
	file <- read.table(paste0(dir, "DepthPCs.vect"), header=F) 
	groups <- read.table(paste0(dir, "Sample.cohort"), header=F)
	uniq.groups <- unique(groups[,2])
	nb.groups <- length(uniq.groups)

	buffer <- 0.01
	xmin <- min(file[,3]) - buffer
	xmax <- max(file[,3]) + buffer
	ymin <- min(file[,4]) - buffer
	ymax <- max(file[,4]) + buffer
	cohort.list<-c("Levine","Davina","Hardcastle","IoO","IoN","IoOFFS","IoONov2013","IoOPanos","Kelsell","LambiaseSD",
	"Lambiase","LayalKC","Manchester","Nejentsev","PrionUnit","Prionb2","Shamima","Sisodiya","Syrris","Vulliamy","WebsterURMD")
	png(paste0(dir, "/plots/DepthPCA.png") ) 
	for(i in 1:nb.groups)
#	for(i in 72:nb.groups)
	{
		hit <- which(groups[,2] == uniq.groups[i])
		if(i==1)
		{
		plot(file[hit,3], file[hit,4], 
			xlab = "PC1", ylab = "PC2", 
			xlim=c(xmin, xmax), 
			ylim=c(ymin, ymax), 
			main = paste("DepthPCA", date()) ,
			col=i
			) 
		} else
		{
			points(file[hit,3], file[hit,4], col=i)
		}
	}
	dev.off()
	png(paste0(dir, "/plots/DepthPCA.var.png") ) 
	file <- read.table(paste0(dir,"DepthPCs.values"), header=F)
	file$sd <- file[,1]^.5
	file$var <- file[,1]^2 / sum(file[,1]^2)
	plot(file$var*100, xlab = "Read Depth Principal Components" , ylab = "Proportion Variance Explained (%)", main = "Variance explained by each Principal Component") 
	dev.off() 

	' >> $oFile
$R CMD BATCH --no-save --no-restore $oFile