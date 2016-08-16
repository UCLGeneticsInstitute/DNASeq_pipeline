getArgs <- function() {
  myargs.list <- strsplit(grep("=",gsub("--","",commandArgs()),value=TRUE),"=")
  myargs <- lapply(myargs.list,function(x) x[2] )
  names(myargs) <- lapply(myargs.list,function(x) x[1])
  return (myargs)
}

rootODir<-'/SAN/vyplab/UCLex/mainset_July2016/cian/'
myArgs <- getArgs()

if ('rootODir' %in% names(myArgs))  rootODir <- myArgs[[ "rootODir" ]]
if ('release' %in% names(myArgs))  release <- myArgs[[ "release" ]]
#############################################

base<-read.table(paste0(rootODir,'Phenotypes'),header=F)
base.groups<-read.table(paste0(rootODir,'GroupNames'),header=F)
fam<-read.table(paste0(rootODir,'allChr_snpStats_out.fam'),header=F)

cohort.list<-c('Levine','Hardcastle','IoO','IoN','Kelsell','LambiaseSD','Lambiase','LayalKC','Nejentsev','PrionUnit','Prionb2','Shamima','Sisodiya','Syrris','Vulliamy','WebsterURMD')
write.table(cohort.list,paste0(rootODir,"cohort.list"),col.names=F,row.names=F,quote=F,sep="\t") 
res<-list.files(paste0(rootODir,'KinshipDecomposition'),pattern="\\.res",full.names=T)

res<-res[grep('tech',res)]
for(i in 1:length(cohort.list))
{

	res.file<-res[grep(cohort.list[i],res)]
	if(cohort.list[i]=="Lambiase") res.file<-res.file[-grep('SD',res.file)]
	if(cohort.list[i]=="IoO") res.file<-res[grep('IoO\\.',res)]

	file<-read.table(res.file,header=T,sep="\t")

	if(i==1)
	{
		dat<-data.frame(matrix(nrow=nrow(file),ncol=(length(cohort.list)+2)))
		residuals<-data.frame(matrix(nrow=nrow(file),ncol=(length(cohort.list)+2)))
		dat[,1:2]<-file[,1:2]
		residuals[,1:2]<-file[,1:2]
	}
	#dat[,i+2]<-file$Phenotype
	#residuals[,i+2]<-file$Residual
	dat[dat[,1]%in%file[,1],i+2]<-file$Phenotype[file[,1]%in%dat[,1]]
	residuals[residuals[,1]%in%file[,1],i+2]<-file$Residual[file[,1]%in%residuals[,1]]
}


perm<-data.frame(matrix(nrow=nrow(file),ncol=(length(cohort.list)+2)))
perm[,1:2]<-dat[,1:2]
perm[,3:ncol(perm)]<-0
for(i in 1:length(cohort.list))
{
	nb.cases<-length(grep(1,dat[,i+2]))
	perm[sample(1:nrow(perm),nb.cases),i+2]<-1
}

ancestry<-read.table(paste0(dir,"UCLex_samples_ancestry"),header=T) 
caucasian<-ancestry[ancestry$Caucasian,1]
dat[!dat[,1]%in%caucasian,3:ncol(dat)]<-NA
perm[!perm[,1]%in%caucasian,3:ncol(perm)]<-NA
residuals[!residuals[,1]%in%caucasian,3:ncol(residuals)]<-NA

#dat<-data.frame(base[,1:2],
exit

write.table(dat,paste0(rootODir,'Clean_pheno_subset'),col.names=F,row.names=F,quote=F,sep="\t")
write.table(perm,paste0(rootODir,'Clean_pheno_subset_permuted'),col.names=F,row.names=F,quote=F,sep="\t")
write.table(residuals,paste0(rootODir,'Clean_pheno_subset_tk_depth_residuals'),col.names=F,row.names=F,quote=F,sep="\t")

dat[dat==1]<-2
dat[dat==0]<-1
dat[is.na(dat)]<- '-9'

write.table(dat,paste0(rootODir,'Clean_pheno_subset_plink'),col.names=F,row.names=F,quote=F,sep="\t")
