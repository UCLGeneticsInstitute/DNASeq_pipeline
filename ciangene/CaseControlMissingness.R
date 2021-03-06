getArgs <- function() {
  myargs.list <- strsplit(grep("=",gsub("--","",commandArgs()),value=TRUE),"=")
  myargs <- lapply(myargs.list,function(x) x[2] )
  names(myargs) <- lapply(myargs.list,function(x) x[1])
  return (myargs)
}


myArgs <- getArgs()

if ('rootODir' %in% names(myArgs))  rootODir <- myArgs[[ "rootODir" ]]
if ('release' %in% names(myArgs))  release <- myArgs[[ "release" ]]

######################
library(parallel) 

plink<-'/home/sejjcmu/bin/plink/plink --noweb --allow-no-sex --bfile' 
outDir<-paste0(rootODir,"CaseControlMissingness/") 
if(!file.exists(outDir))dir.create(outDir) 

groups<-read.table(paste0(rootODir,"GroupNames"),header=F) 
nb.groups<-nrow(groups) 

for(i in 1:nrow(groups))
{
	oFile<- paste0(outDir,groups[i,1],".sh")  
	write.table(paste0("release=",release),oFile,col.names=F,row.names=F,quote=F,sep="\t") 
	write.table(paste0("rootODir=",rootODir),oFile,col.names=F,row.names=F,quote=F,sep="\t",append=T) 
	write.table(paste0("mpheno=",i),oFile,col.names=F,row.names=F,quote=F,sep="\t",append=T) 
	write.table(paste0("batch=",groups[i,1]),oFile,col.names=F,row.names=F,quote=F,sep="\t",append=T) 
	system(paste('cat /SAN/vyplab/UCLex/scripts/DNASeq_pipeline/ciangene/CaseControlMissingnessTemplate >>', oFile)) 

}
files<-list.files(outDir,pattern="sh",full.names=T)
names<-gsub(basename(files),pattern="\\..*",replacement="")
cohort.list<-c('Levine','Hardcastle','IoO','IoN','Kelsell','LambiaseSD','Lambiase','LayalKC','Nejentsev','PrionUnit','Prionb2','Shamima','Sisodiya','Syrris','Vulliamy','WebsterURMD')
files<-files[names%in%cohort.list]
mclapply(files,function(x)system(paste("sh",x)),mc.cores=2) 
