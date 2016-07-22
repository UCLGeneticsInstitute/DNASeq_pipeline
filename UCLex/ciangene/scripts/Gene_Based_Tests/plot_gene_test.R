getArgs <- function() {
  myargs.list <- strsplit(grep("=",gsub("--","",commandArgs()),value=TRUE),"=")
  myargs <- lapply(myargs.list,function(x) x[2] )
  names(myargs) <- lapply(myargs.list,function(x) x[1])
  return (myargs)
}

release <- 'February2015'
rootODir<-'/scratch2/vyp-scratch2/cian'
myArgs <- getArgs()

if ('rootODir' %in% names(myArgs))  rootODir <- myArgs[[ "rootODir" ]]
if ('release' %in% names(myArgs))  release <- myArgs[[ "release" ]]

######################

library(snpStats)

oDir <- paste0("/scratch2/vyp-scratch2/cian/UCLex_", release, "/") 

files <- list.files(paste0(oDir, "LDAK_gene_tests_all_phenos/"), pattern = "with_kin", full.names=T) 
#batches <- gsub(files, pattern = ".*kin_", replacement = "")
batches <- basename(files)
names <- gsub(batches, pattern = "with_kin_", replacement = "") 
cohorts <- gsub(names, pattern ="_.*", replacement = "") 	

groups <- read.table(paste0(oDir, "cohort.summary"), header=T) 
groups$GenomicInflationTech <- 0 
groups$GenomicInflationPerm <- 0 
groups$GenomicInflationBase <- 0 

noKins <- paste0(oDir, "LDAK_gene_tests_all_phenos/no_kin", gsub(batches,pattern="with_kin", replacement=''),  "/regressALL")
res <- paste0(oDir, "LDAK_gene_tests_all_phenos/no_kin", gsub(batches,pattern="with_kin", replacement=''),  "_res/regressALL")

var <- list.files(paste0(oDir, "KinshipDecomposition") , pattern = "tech.progress", full.names=T) 


oFile <- paste0( "UCLex_", release, "_gene_tests.pdf") 
pdf(oFile)
par(mfrow=c(2,2), cex.main = 0.6) 

	for(i in 1:length(files)) 	
	{
		target <- grep(paste0(cohorts[i],"_"), var)
		if(length(target)>0)
		{
			cohort.res <- read.table(var[grep(paste0(cohorts[i],"_"), var)  ], header=T, sep="\t") 
			var.explained <- paste0("(", round(cohort.res[nrow(cohort.res),2], 3 ) *100, "%)")	
		}
		if(file.exists(noKins[i]))	
		{
			file <- read.table(noKins[i] , header=T,stringsAsFactors=F)		
			GI <- qq.chisq(-2*log(as.numeric(file$LRT_P_Perm)), df=2, x.max=30, main = paste(names[i], "Base", sum(file$Gene_Weight)) ,  pvals=T)
			legend("topleft", , signif(GI[3][[1]]),2) 
			out <- gsub(dirname(noKins[i]), pattern = ".*phenos/", replacement = "")
			message(paste(out, "base GI is", GI[3][[1]] ) ) 

			hit <- grep(batches[i], groups$Cohort ) 		
			groups$GenomicInflationBase[hit] <- GI[3][[1]]
		}
		permFile <- paste0(files[i], "/regress1")
		if(file.exists(permFile))
		{
			file <- read.table(permFile , header=T,stringsAsFactors=F)	
			pvals=pchisq(file$Perm_1,1,lower=F)*0.5
			GI <- qq.chisq(-2*log(as.numeric(pvals)), df=2, x.max=30, main = paste(names[i], "Perm", sum(file$Gene_Weight)) ,  pvals=T)
			legend("topleft", , signif(GI[3][[1]]),2) 
			out <- gsub(dirname(permFile), pattern = ".*phenos//", replacement = "")
			message(paste(out, "permuted GI is", GI[3][[1]] ) ) 

			hit <- grep(batches[i], groups$Cohort ) 		
			groups$GenomicInflationPerm[hit] <- GI[3][[1]]
		}
		techFile <- paste0(files[i], "/regressALL")
		if(file.exists(techFile)) 
		{	
			file <- read.table(techFile , header=T,stringsAsFactors=F)		
			GI <- qq.chisq(-2*log(as.numeric(file$LRT_P_Perm)), df=2, x.max=30, main = paste(names[i], "Tech", sum(file$Gene_Weight), var.explained) ,  pvals=T)
			legend("topleft", , signif(GI[3][[1]]),2) 
			out <- gsub(dirname(techFile), pattern = ".*phenos//", replacement = "")
			message(paste(out, "GI is", GI[3][[1]] ) ) 

			hit <- grep(batches[i], groups$Cohort ) 		
			groups$GenomicInflationTech[hit] <- GI[3][[1]]
		}
		if(file.exists(res[i])) 
		{	
			file <- read.table(res[i], header=T,stringsAsFactors=F) 		
			GI <- qq.chisq(-2*log(as.numeric(file$LRT_P_Perm)), df=2, x.max=30, main = paste(names[i], "TechResiduals", sum(file$Gene_Weight), var.explained) ,  pvals=T)
			legend("topleft", , signif(GI[3][[1]]),2) 
			out <- gsub(dirname(techFile), pattern = ".*phenos//", replacement = "")
			message(paste(out, "GI is", GI[3][[1]] ) ) 

			hit <- grep(batches[i], groups$Cohort ) 		
			groups$GenomicInflationTech[hit] <- GI[3][[1]]
		}
	}	
dev.off()



