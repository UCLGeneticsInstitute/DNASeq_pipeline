plink<-'/share/apps/genomics/plink-1.07-x86_64/plink --noweb --allow-no-sex --bfile'

removeConflictingControls<-function(basePheno,remove,cases,oDir)
{
	# basePheno is initial pheno file, eg the one from make_phenotype_file.R

	groups<-read.table(paste0(oDir,"cohort.summary"),header=T)
	group<-gsub(cases[1],pattern="_.*",replacement='')
	
	case.rows<-basePheno[basePheno[,1]%in% cases,]
	case.col<-which(basePheno[basePheno[,1]%in% cases,][1,]==2)
	if(length(case.col)>1)case.col<-which(groups$Cohort %in% group) + 2

	case.coldata<-basePheno[,case.col]
	for(i in 1:length(remove))
	{
		if(i==1)filtPheno<-basePheno
		hit<-grep(remove[i],basePheno[,1])
		if(length(hit)==0){ warning(paste("No samples called" ,remove[i],"found"))}else{
			filtPheno[hit,case.col]<-NA
		}
	}
return(filtPheno)
}

makeExternalControls<-function(pheno,cases,data,oBase,percent=10)
{
	# pheno is phenotype file, preferably with any conflicting controls removed. 
	# list of case names
	# percent is what percent of controls i want to use as my external control set. 
	# data is stem for bam file 
	out<- paste0(oBase,"_ex_ctrls")
	
	case.rows<-pheno[grep(cases, pheno[,1] ),]
	case.col<-colnames(pheno)%in%cases 

	duds<-which(is.na(pheno[,case.col]))
	ctrls<- pheno[!pheno[,1] %in% case.rows[,1] ,1]
	if(length(duds)>0)ctrls<-ctrls[!ctrls%in%pheno[duds,1]] ; print(paste(cases, 'has', length(duds), 'missing samples') ) 
	ctrls<- ctrls[!is.na(ctrls)]
	ctrls<-ctrls[-grep("One",ctrls)]
	
	if(file.exists(out))
	{
		print("External Controls already exist, so just fixing phenotype file directly") 
		ex.ctrls<-read.table(out,header=T,sep='\t') 
		pheno[ pheno[,1]%in%ex.ctrls , case.col] <- NA
	} else
	{
		nb.ex.ctrl<-round(length(ctrls)/percent) 
		ex.ctrls<-ctrls[sample(1:length(ctrls),nb.ex.ctrl)]	
		pheno[ pheno[,1]%in%ex.ctrls , case.col] <- NA

		write.table(data.frame(ex.ctrls,ex.ctrls),out,col.names=F,row.names=F,quote=F)
		run<-paste( plink,data, '--freq --out', out ,'--keep', out) 
		#system(run)
		run2<-paste( plink,data, '--hardy --out',out ,'--keep', out ) 
		#system(run)
		run3<-paste( plink,data, '--missing --out',out ,'--keep', out ) 
		#system(run)
		runs<-c(run,run2,run3)
		mclapply(runs,function(x)system(x),mc.cores=3)
	}
	out<-paste0(oBase,'_CC_controls')
	if(!file.exists(out))
	{
		ctrls<-ctrls[!ctrls %in% ex.ctrls]
		write.table(data.frame(ctrls,ctrls),out,col.names=F,row.names=F,quote=F)
		run<-paste( plink,data, '--freq --out', out ,'--keep', out) 
		#system(run)
		run2<-paste( plink,data, '--hardy --out',out ,'--keep', out ) 
		#system(run)
		run3<-paste( plink,data, '--missing --out',out ,'--keep', out ) 
		#system(run)
		runs<-c(run,run2,run3)
		mclapply(runs,function(x)system(x),mc.cores=3)
	}
	return(pheno)
}


