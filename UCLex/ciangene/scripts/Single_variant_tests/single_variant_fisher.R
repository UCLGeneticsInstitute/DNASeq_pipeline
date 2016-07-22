getArgs <- function() {
  myargs.list <- strsplit(grep("=",gsub("--","",commandArgs()),value=TRUE),"=")
  myargs <- lapply(myargs.list,function(x) x[2] )
  names(myargs) <- lapply(myargs.list,function(x) x[1])
  return (myargs)
}

release <- 'February2015'

myArgs <- getArgs()

if ('rootODir' %in% names(myArgs))  rootODir <- myArgs[[ "rootODir" ]]
if ('release' %in% names(myArgs))  release <- myArgs[[ "release" ]]

######################

library(MultiPhen)

bDir <- paste0("/scratch2/vyp-scratch2/cian/UCLex_", release ) 

files <- list.files(bDir, pattern = "bim", full.names=T)
iFile <- gsub(files, pattern = "\\.bim", replacement = "")
iFile <- iFile[grep("rol/UCL", iFile) ]

pheno <- read.table(pheno, header=F, sep="\t")
pheno <- subset(pheno, pheno[,3] != '-9') 

plink <- read.plink(iFile[fil])
bim <- read.table( paste0(iFile[fil], '.bim'), header=F, sep="\t")
colnames(plink) <- bim[,2]

colNamesDat  <- c("SNP", "FisherPvalue", "nb.mutations.cases", "nb.mutations.ctrls", "nb.cases", "nb.ctrls", "case.maf", "ctrl.maf", 
				"Case.Hwe", "Ctrl.Hwe", "CC.Hwe", "CC.maf", "nb.NAs.cases", "nb.NAs.ctrls", "nb.Homs.cases", "nb.Homs.ctrls"
	 				,"case.call.rate", "ctrl.call.rate"
	 			) 
dat <- data.frame(matrix(nrow = ncol(plink), ncol = length(colNamesDat) ) )
colnames(dat) <- colNamesDat
dat[,1] <- colnames(plink) 

library(snpStats) ## this is down here because its "read.plink" function works differently to MultiPhens

for(i in 1:ncol(plink))
{
	case.calls <- t(data.frame(plink[ pheno[,phenoType+2 ] == 1 , i])) 
	rownames(case.calls) <- colnames(plink)[i]
	ctrl.calls <- t(data.frame(plink[ pheno[,phenoType+2] == 0 , i]) )
	rownames(ctrl.calls) <- colnames(plink)[i]

	number_mutations_cases <- sum( case.calls , na.rm=T )
	number_mutations_ctrls <- sum( ctrl.calls , na.rm=T ) 

	number_Homs_cases <- length(which(unlist(case.calls) == 2))
	number_Homs_ctrls <- length(which(unlist(ctrl.calls) == 2))
	
	nb.nas.cases <- length(which(is.na(case.calls)))
	nb.nas.ctrls <- length(which(is.na(ctrl.calls)))

	nb.cases <-  length(which(!is.na( case.calls )) ) 
	nb.ctrls <-  length(which(!is.na( ctrl.calls )) ) 
		
	mean_number_case_chromosomes <- nb.cases * 2
	mean_number_ctrl_chromosomes <- nb.ctrls * 2

	if (!is.na(number_mutations_cases)  & !is.na(number_mutations_ctrls)  )
	{
	if (nb.cases > 0 & nb.ctrls > 0)
	{
		fishertest <-  fisher.test((matrix(c(number_mutations_cases, mean_number_case_chromosomes
		                         - number_mutations_cases, number_mutations_ctrls, mean_number_ctrl_chromosomes - number_mutations_ctrls),
		                       nrow = 2, ncol = 2)))


		dat$FisherPvalue[i] 		<- fishertest$p.value 
		dat$nb.mutations.cases[i] 	<- number_mutations_cases
		dat$nb.mutations.ctrls[i] 	<- number_mutations_ctrls
		dat$nb.cases[i]				<- nb.cases 
		dat$nb.ctrls[i] 			<- nb.ctrls 
		dat$nb.NAs.cases[i] 		<- nb.nas.cases 
		dat$nb.NAs.ctrls[i] 		<- nb.nas.ctrls
		dat$nb.Homs.cases[i]		<- number_Homs_cases
		dat$nb.Homs.ctrls[i]		<- number_Homs_ctrls	
	}
	}
}


chr <- gsub(iFile[fil], pattern = ".*_", replacement = "") 

release <- 'October2014'
dir <- paste0("/cluster/project8/vyp/exome_sequencing_multisamples/mainset/GATK/mainset_", release , "/mainset_", release, "_by_chr")
files <- list.files(dir, pattern ="_snpStats.RData", full.names=T)
if(as.numeric(chr) == 23) chr <- 'X'

current.file <- files[grep(paste0("chr", chr, "_") , files) ]
load(current.file)

#mat.snp.cols <- t(data.frame(strsplit(as.character(colnames(matrix.calls.snpStats)) , "_")))
#mat.snp.cols <- paste0(mat.snp.cols [,1], "_", mat.snp.cols [,2] ) 

matrix.calls.snpStats <- matrix.calls.snpStats[ , colnames(matrix.calls.snpStats)  %in% dat[,1] ]

sample.qc <- row.summary(matrix.calls.snpStats) 

case.rows <- rownames(sample.qc) %in% colnames(case.calls)  
calls.cases <- matrix.calls.snpStats[unlist(case.rows) ,]
case.summary <- col.summary(calls.cases)
case.maf <- case.summary$MAF
case.hwe <- case.summary$z.HWE
case.call.rate <- case.summary$Call.rate

ctrl.rows <- rownames(sample.qc) %in% colnames(ctrl.calls)  
calls.ctrls <- matrix.calls.snpStats[unlist(ctrl.rows) ,]
ctrl.summary <- col.summary(calls.ctrls)
ctrl.maf <- ctrl.summary$MAF
ctrl.hwe <- ctrl.summary$z.HWE
ctrl.call.rate <- ctrl.summary$Call.rate

CC.calls <- matrix.calls.snpStats[c(unlist(which(case.rows)), unlist(which(ctrl.rows))), ]
CC.summary <- col.summary(CC.calls)
CC.hwe <- CC.summary$z.HWE
CC.maf <- CC.summary$MAF

dat <- dat[ dat[,1] %in% colnames(matrix.calls.snpStats) , ]
dat$case.maf <- case.maf
dat$ctrl.maf <- ctrl.maf
dat$Case.Hwe <- case.hwe
dat$Case.Hwe <- ctrl.hwe
dat$CC.Hwe <- CC.hwe
dat$CC.maf <- CC.maf

dat$case.call.rate <- case.call.rate
dat$ctrl.call.rate <- ctrl.call.rate


oDir <- "/cluster/project8/vyp/cian/data/UCLex/UCLex_October2014/All_phenotypes/Fisher_single_variant/"

chr <- as.numeric(gsub(rownames(case.calls)[1], pattern ="_.*", replacement = "") )
CC.out <- paste0(oDir, oFile, "_Chr_", chr) 

write.table(dat, CC.out , col.names=T, row.names=F, quote= F, sep="\t", append =F)
