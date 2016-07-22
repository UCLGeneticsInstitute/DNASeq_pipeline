comp <- "mbsp"
if(comp == "mbp") options(width=170)

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
#############################################

oDir <- paste0(rootODir, "/UCLex_", release, "/")

fam <- read.table(paste0(oDir, "UCLex_",release,".fam"), header=F) 
groups <- gsub(fam[,1], pattern = "_.*",replacement = "")
groups.unique <- unique(groups)
# Tring to group samples by cohort correctly. Default method is to use string before first underscore in their name, but that doesn't work for all samples, 
# so to try group correctly (eg to prevent ALevine from being treated as a separate phenotype from Levine, use fixPhenoGroupings )
use.fixPhenoGroupings <- TRUE

if(use.fixPhenoGroupings)
{
	cohorts.to.fix <- c("Levine", "B2", "BC", "UCL", "UCLG", "Syrris", "gosgeneBGI")
	
	nb.original.groups <- length(groups.unique)

	fixPhenoGroupings <- function(cohorts, inGroups)
	{

		for(i in 1:length(cohorts.to.fix))
		{
			hit <- grep(cohorts.to.fix[i], inGroups)
			inGroups[hit] <- cohorts.to.fix[i]
		}

		groups.unique <- unique(inGroups)
		nb.new.groups <- length(groups.unique)

		out <- paste( nb.original.groups - nb.new.groups , "samples have been merged into other cohorts" )
		message(out)

	return(inGroups)
	} # fixPhenoGroupings 

	groups <- fixPhenoGroupings(cohorts.to.fix, groups)
	groups.unique <- unique(groups)


} ### use.fixPhenoGroupings 

nb.groups <- length(unique(groups))

## now make phenotype file 
pheno <- data.frame(matrix(nrow = nrow(fam), ncol = nb.groups+2))
colnames(pheno) <- c("SampleID1", "SampleID2", groups.unique)
pheno[,1:2] <- fam[,1]

for(i in 1:nb.groups)
{
	pheno[, i + 2] <- 1 
	hits <- grep(groups.unique[i], groups )
#	message(length(hits) 	) 
	pheno[hits,i+2] <- 2
}



## remove pheno for extCtrls
#extCtrls <- read.table(paste0(oDir, "ext_ctrl_samples"), header=F) ## made in first.step.R
#ex.ctrl.pheno <- pheno[,1] %in% unlist(extCtrls) 
#pheno[ex.ctrl.pheno,3:ncol(pheno)] <- 'NA'

## remove Lambiase family
family <- c("UCLG569", "UCLG567", "LambiaseSD_UCLG594", "UCLG568", "UCLG570", "UCLG571", "UCLG572" )
for(i in 1:length(family))
{
	hit <- grep(family[i], pheno[,1])
	pheno[hit,3:ncol(pheno)] <- 'NA' 
} 

## summarise case cohort breakdowns
cohort.summary <- data.frame(do.call(rbind, lapply(pheno[,3:ncol(pheno)], table) )) 
cohort.summary$Cohort <- rownames(cohort.summary) 
colnames(cohort.summary) <- c("Nb.Ctrls", "Nb.cases", "Nb.ext.Ctrls", "Cohort") 
head(cohort.summary)

for(i in 1:nb.groups)
{
	cohort.breakdown<- data.frame(table(as.character(pheno[i,3:ncol(pheno)]) ) )
	nb.ctrls<-grep(1,cohort.breakdown[,1])
	nb.cases<-grep(2,cohort.breakdown[,1])
	nb.nas<-grep('NA',cohort.breakdown[,1])

	if(length(nb.ctrls)>0)cohort.summary$Nb.Ctrls[i]<-cohort.breakdown[nb.ctrls,2] 
	if(length(nb.cases)>0)cohort.summary$Nb.cases[i]<-cohort.breakdown[nb.cases,2] 
	if(length(nb.nas)>0)cohort.summary$Nb.ext.Ctrls[i]<-length(grep(pheno[i,1],pheno[,1])) 
}


write.table(groups.unique, paste0(oDir, "GroupNames"), col.names=F, row.names=F, quote=F, sep="\t")
write.table(pheno, file = paste0(oDir, "Phenotypes"), col.names=F, row.names=F, quote=F, sep="\t") 
write.table(data.frame(pheno[,1], groups ), paste0(oDir, "Sample.cohort"),  col.names=F, row.names=F, quote=F, sep="\t") 
write.table(cohort.summary, file = paste0(oDir, "cohort.summary"), col.names=T, row.names=F, quote=F, sep="\t") 


### fastLMM wants missing as -9, so use separate pheno file. 
oDir <- paste0(rootODir, "/UCLex_", release, "/")
inPheno<-paste0(oDir, 'Phenotypes')
pheno <- read.table(inPheno, header=F)
pheno[is.na(pheno)] <- '-9'
write.table(pheno,paste0(inPheno,"_fastlmm"), col.names=F, row.names=F, quote=F, sep="\t")

###### make permuted phenotype file for fastLMM/plink 
for(i in 1:nb.groups)
{
	nb.cases <- length(which(pheno[,i+2]==2)) 
	pheno[,i+2] <- 1 
	pheno[sample(1:nrow(pheno),nb.cases),i+2] <- 2
}
write.table(pheno,paste0(inPheno,"_fastlmm_permuted"), col.names=F, row.names=F, quote=F, sep="\t")



