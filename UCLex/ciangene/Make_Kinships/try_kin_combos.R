kins<-list.files(pattern='details')
kins<-gsub(kins,pattern="\\.grm.details",replacement="")
tk<-kins[grep("TK",kins)]
rd<-kins[grep("Depth",kins)]
combos<-expand.grid(tk,rd)
runSh='sh /cluster/project8/vyp/cian/scripts/bash/runBashCluster.sh'
dir<-'/scratch2/vyp-scratch2/cian/UCLex_July2015/'
ldak<-'/cluster/project8/vyp/cian/support/ldak/ldak'
oDir<-paste0(dir,'Combo/')
if(!file.exists(oDir))dir.create(oDir)
phenotype<-paste0(dir,'Clean_pheno_subset')
mpheno<-7

prep<-FALSE
if(prep)
{
	kinDir<-paste0(getwd(),'/')
	for(i in 1:nrow(combos))
	{
		kinshipList<-paste0(oDir,'kinships_',i) 
		tk.kin<-paste0(kinDir,combos[i,1])
		rd.kin<-paste0(kinDir,combos[i,2])
		write.table(data.frame(tk.kin,rd.kin),kinshipList,col.names=F,row.names=F,quote=F,sep="\t") 
		oFile<-paste0(oDir,'run',i,'.sh')
		for(pheno in 1:16)
		{
			run<-paste(ldak,'--reml',paste0(oDir,'kins',i,'_',pheno),'--mgrm',kinshipList,'--pheno',phenotype,'--mpheno',pheno,'--time-save YES') 

			write.table(run,oFile,col.names=F,row.names=F,quote=F,sep="\t",append=T)
		}
		system(paste(runSh,oFile))
	}
}


test<-TRUE
if(test)
{
	files<-list.files(oDir,pattern='indi.res',full.names=T)
	fastlmm='/share/apps/bahler/FaSTLMM.207/Bin/Linux_MKL/fastlmmc'
	data=paste0(dir,'allChr_snpStats_out')
	runSh='sh /cluster/project8/vyp/cian/scripts/bash/runBashCluster_large.sh'
	batches<-unique(gsub(basename(files),pattern='_.*',replacement=''))
	for(i in 1:length(batches))
	{
		hit<-grep(paste0(batches[i],'_'),files)
		if(length(hit)>0)
		{
		oFile<-paste0(oDir,'test',i,'.sh')
		file.remove(oFile)
		for(hi in 1:length(hit))
		{
			file<-read.table(files[hi],header=T,sep="\t")
			pheno<-data.frame(file[,1:2],file$Residual)
			oPheno<-paste0(oDir,batches[i],'pheno',hi)
			write.table(pheno,oPheno,col.names=F,row.names=F,quote=F,sep="\t",append=F)
			run<-paste(fastlmm,'-linreg -simLearnType Full -verboseOutput -numJobs 10 -thisJob 1 -missingPhenotype NA -maxThreads 1 -bfile',data,'-pheno',oPheno,'-out',paste0(oDir,batches[i],'test',i ) )  
			write.table(run,oFile,col.names=F,row.names=F,quote=F,sep="\t",append=T)
		}
		system(paste(runSh,oFile))
		}
	}
}

process<-FALSE
if(process)
{
	files<-list.files(oDir,pattern='test',full.names=T)
	files<-files[-grep("sh",files)]

	levine<-read.table(paste0(dir,'FastLMM_Single_Variant_all_phenos/Levine_final'),header=T,sep="\t")
	base<-data.frame(SNP=levine$SNP,Pvalue=as.numeric(as.character(levine$Pvalue) ) ) 

	results<-data.frame(matrix(nrow=length(files),ncol=2))
	results[,1]<-basename(files)
	for(i in 1:length(files))
	{
		file<-read.table(files[i],header=T,sep="\t")
		pvals<-data.frame(SNP=file$SNP,pval=as.numeric(as.character(file$Pvalue) ) ) 
		dat<-merge(pvals,base,by='SNP')
		results[i,2]<-cor(dat$pval,dat$Pvalue)
		if(i==1)best<-results[i,2]
		if(results[i,2]>best)
		{
			best<-results[i,2]
			print(paste('New best is',best))
		}	
	}
	write.table(results,'results',col.names=T,row.names=F,quote=F,sep="\t")
}

final<-FALSE
if(final)
{
	results<-read.table('results',header=T) 
	top<-results[ which(results[,2]==max(results[,2])),] 
	top$Count<-gsub(top[,1],pattern="test",replacement='')
	kins<-paste0(oDir,'/kinships_',top$Count)
	lapply(kins,function(x)system(paste('cat',x,'>>',paste0(oDir,'merged')))) 

	lapply(kins,function(x)system(paste('cp',x,'>>',paste0(oDir,x,'_nice')))) 
}





