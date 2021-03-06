getArgs <- function() {
  myargs.list <- strsplit(grep("=",gsub("--","",commandArgs()),value=TRUE),"=")
  myargs <- lapply(myargs.list,function(x) x[2] )
  names(myargs) <- lapply(myargs.list,function(x) x[1])
  return (myargs)
}

myArgs <- getArgs()

if ('rootODir' %in% names(myArgs))  rootODir <- myArgs[[ "rootODir" ]]
if ('release' %in% names(myArgs))  release <- myArgs[[ "release" ]]
#############################################
library(parallel) 
source("/SAN/vyplab/UCLex/scripts/DNASeq_pipeline/ciangene/MakePhenotypes/make_phenotype_file.R")
source("/SAN/vyplab/UCLex/scripts/DNASeq_pipeline/ciangene/MakePhenotypes/CaseControl_support.R")

groups<-read.table( paste0(rootODir, "cohort.summary"), header=T,sep="\t")

data=paste0(rootODir,"allChr_snpStats_out")
inPheno<-paste0(rootODir,"Phenotypes")
pheno<-read.table(inPheno,header=F,sep="\t") 
colnames(pheno) <- c("SampleID1", "SampleID2", groups[,ncol(groups)])
remove.caucasians<-TRUE
if(remove.caucasians)
{
	ancestry<-read.table(paste0(rootODir,"UCLex_samples_ancestry"),header=T,sep="\t")
	remove<-subset(ancestry,ancestry$Caucasian==FALSE)
	pheno[pheno[,1]%in%remove[,1],3:ncol(pheno)]<-NA
}

syrris<-removeConflictingControls(pheno,remove=c("Lambiase") ,cases=pheno[grep("Syrris",pheno[,1]),1],oDir=rootODir) 
lambiase<-removeConflictingControls(syrris,remove=c("Syrris"),cases=pheno[grep("Lambiase_",pheno[,1]),1] ,oDir=rootODir ) 
lambiaseSD<-removeConflictingControls(lambiase,remove=c("Syrris","Lambiase_"),cases=pheno[grep("LambiaseSD",pheno[,1]),1] ,oDir=rootODir ) 
##lambiaseSD[c(grep("yrris",lambiaseSD[,1]),grep("ambiase",lambiaseSD[,1]) ),c(1,94,95,105)] ## check, syrris should be NA in cols 92 93, lambiase  and syrris NA in 92 and 94 and SD only na in 104. 
rows<-c(grep("ambiase",pheno[,1]),grep("yrris",pheno[,1]))

runSh='sh /cluster/project8/vyp/cian/scripts/bash/runBashCluster.sh'

cohort.list<-c('Levine','Davina','Hardcastle','IoO','IoN','IoOFFS','IoONov2013','IoOPanos','Kelsell','LambiaseSD',
'Lambiase','LayalKC','Manchester','Nejentsev','PrionUnit','Prionb2','Shamima','Sisodiya','Syrris','Vulliamy','WebsterURMD')

if(!file.exists(paste0(rootODir,"/External_Control_data/") ))dir.create(paste0(rootODir,"External_Control_data/") ) 
pheno<-lambiaseSD  ## or last modified pheno file
data.frame( pheno[rows,1], pheno$Lambiase[rows],pheno$LambiaseSD[rows],pheno$Syrris[rows] ) 

plink<-'/home/sejjcmu/bin/plink/plink --noweb --allow-no-sex --bfile' 
for(i in 1:length(cohort.list))
{
	hit<-grep(cohort.list[i],groups$Cohort)
	if(cohort.list[i] == "IoO")	hit<-grep('IoO$',groups$Cohort)
	if(cohort.list[i] == "Lambiase") 	hit<-grep('Lambiase$',groups$Cohort)
	pheno<-makeExternalControls(pheno,cases=groups$Cohort[hit],data=data,oBase=paste0(rootODir,"External_Control_data/",groups$Cohort[hit]) ) 
	case.samples<-pheno[grep(cohort.list[i],pheno[,1]),1]
	case.list<-paste0(rootODir,"External_Control_data/",groups$Cohort[hit],"_case_samples")
	write.table(data.frame(case.samples,case.samples),  case.list ,col.names=F,row.names=F,quote=F,sep="\t") 
	run<-paste(plink, data, '--missing --out', paste0(rootODir,"External_Control_data/",groups$Cohort[hit],"_case_qc") , '--keep',case.list) #
#	system(run) 
	run2<-paste(plink, data, '--freq --out', paste0(rootODir,"External_Control_data/",groups$Cohort[hit],"_case_maf")  , '--keep',case.list) #
	write.table(run, paste0(rootODir,"External_Control_data/",groups$Cohort[hit],".sh") ,col.names=F,row.names=F,quote=F,sep="\t") 
	write.table(run2, paste0(rootODir,"External_Control_data/",groups$Cohort[hit],".sh") ,col.names=F,row.names=F,quote=F,sep="\t",append=T) 
##	system(paste(runSh,paste0(oDir,"/External_Control_data/",groups$Cohort[hit],".sh") ) ) 
#	system(run2) 
}

files<-list.files( paste0(rootODir,"External_Control_data/"),pattern="sh",full.names=T)
test<-list.files(paste0(rootODir,"External_Control_data/"),pattern="case_maf")
if(length(test)==0)mclapply(files,function(x)system(paste("sh",x)),mc.cores=4) else print("Cases already processed") 


write.table(pheno,inPheno, col.names=F, row.names=F, quote=F, sep="\t")
## fastLMM wants missing as -9, so use separate pheno file. 
pheno[is.na(pheno)] <- '-9'
write.table(pheno,paste0(inPheno,"_fastlmm"), col.names=F, row.names=F, quote=F, sep="\t")

perm<-pheno
###### make permuted phenotype file for fastLMM/plink 
for(i in 1:nb.groups)
{
	nb.cases <- length(which(perm[,i+2]==2)) 
	perm[,i+2] <- 1 
	perm[sample(1:nrow(perm),nb.cases),i+2] <- 2
}
write.table(perm,paste0(inPheno,"_fastlmm_permuted"), col.names=F, row.names=F, quote=F, sep="\t")

#source("/SAN/vyplab/UCLex/scripts/DNASeq_pipeline/ciangene/MakePhenotypes/MakeGoodPhenotypeFile.R") 


