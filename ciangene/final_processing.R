getArgs <- function() {
  myargs.list <- strsplit(grep("=",gsub("--","",commandArgs()),value=TRUE),"=")
  myargs <- lapply(myargs.list,function(x) x[2] )
  names(myargs) <- lapply(myargs.list,function(x) x[1])
  return (myargs)
}

release <- 'July2015'
rootODir<-'/scratch2/vyp-scratch2/cian'

myArgs <- getArgs()

if ('rootODir' %in% names(myArgs))  rootODir <- myArgs[[ "rootODir" ]]
if ('release' %in% names(myArgs))  release <- myArgs[[ "release" ]]

######################
library(data.table)
library(xtable)
oDir <- paste0(rootODir, "/UCLex_", release, "/")

singles<-list.files(paste0(oDir,'FastLMM_Single_Variant_all_phenos/'),pattern="csv",full.names=T)
genes<-list.files(paste0(oDir,'LDAK_gene_tests_all_phenos_good/'),full.names=T,recursive=T,pattern="regressALL") 
genes<-genes[grep("0000001_.00000001",genes)]
base<-genes[grep("base",genes)]
base.list<-base[-grep("perm",base)]
fixed<-genes[grep("depthkin",genes)]
names<-gsub(gsub(basename(dirname(base.list)) ,pattern='base_',replacement=''),pattern="_.*",replacement="" )
prep<-TRUE
if(prep)
{
outDir<-paste0(oDir,"CaseControlResults/")
if(!file.exists(outDir)) dir.create(outDir)

for(i in 1:length(fixed))
{
	print(names[i])
	file<-read.table(fixed[i],header=T,sep=" ")
	dat<-data.frame(matrix(nrow=nrow(file),ncol=3)) 
	colnames(dat)<-c("Gene","UncorrectedPvalue","CorrectedPvalue")
	dat$Gene<-file$Gene_name
	dat$CorrectedPvalue<-file$LRT_P_Perm
	base<-read.table(base.list[i],header=T,sep=" ")
	dat$UncorrectedPvalue<-base$LRT_P_Perm
	write.table(dat,paste0(outDir,names[i],'_gene_list.csv'),col.names=T,row.names=F,quote=T,sep=",")
}


lapply(singles,function(x) system(paste('cp',x,outDir)))

geneDir<-'/cluster/project8/vyp/cian/data/Support/CandidateGenes/'
geneInfo<-read.table(paste0(geneDir,'cohort_phenotype.txt'),header=F)

files<-list.files(outDir,full.names=T,pattern='csv')
cohorts<-gsub(basename(files),pattern='_.*',replacement='')
for(i in 1:length(files))
{
	file<-fread(files[i],header=T)
	hit<-which(geneInfo[,1]%in%cohorts[i])#grep(cohorts[i],geneInfo[,1])
	genelist<-paste0(geneDir,'candidate_genes_',geneInfo[hit,2])
	file$Candidate<-FALSE
	if(file.exists(genelist))
	{
		genes<-unlist(read.table(genelist,header=F,sep="\t"))
		if(length(grep('external_gene_name',colnames(file)))==1) file$Candidate[file$external_gene_name%in%genes]<-TRUE
		if(length(grep('external_gene_name',colnames(file)))==0) file$Candidate[file$Gene%in%genes]<-TRUE
		print(head(file)) 
		file<-file[order(-file$Candidate),]
		write.table(file,files[i],col.names=T,row.names=F,quote=T,sep=',')
	}
	write.table(file,files[i],col.names=T,row.names=F,quote=T,sep=',')
}

}

dir<-'/scratch2/vyp-scratch2/cian/UCLex_July2015/CaseControlResults/'
snps<-list.files(dir,pattern='single',full.names=T) 
snps<-snps[grep('csv',snps)]
names<-gsub(basename(snps),pattern="_.*",replacement='')
info<-read.table('/cluster/project8/vyp/cian/data/Support/CandidateGenes/cohort_phenotype.txt',header=F)
pval<-0.000001
for(i in 1:length(snps))
{
	phenotype<-info[info[,1,]%in%names[i],2]
	file<-read.table(snps[i],header=T,sep=',')
	dat<-subset(file,file$Candidate & ( file$TechKinPvalue<=pval|file$FisherPvalue<=pval ) ) 

	print(paste(data.frame(basename(snps[i]),nrow(dat),phenotype)) ) 
	if(nrow(dat)>0)
	{
		if(nrow(dat)>20)dat<-dat[1:20,]
		rownames(dat)<-1:nrow(dat)
		mafs<-dat[,grep('maf',colnames(dat))]
		case.maf<-mafs[,-grep('ctrl',colnames(mafs))]
		ctrl.maf<-mafs[,grep('ctrl',colnames(mafs))]

		filt<-data.frame(SNP=dat$SNP,Func=dat$ExonicFunc, Gene=dat$external_gene_name,tkrdPvalue=dat$TechKinPvalue,
			Fisher=dat$FisherPvalue, CaseMAF=case.maf,CtrlMAF=ctrl.maf, ESP6500maf=dat$ESP6500si_ALL
		)
	filt.xtable<-xtable(filt,caption=paste(phenotype,"Single Variant Results in Candidate Genes") ,digits=2, display = c(rep("s",4),rep("E",5)))
	print.xtable(filt.xtable, type="latex",file=paste0(dir,names[i],"_snp.tex"),scalebox=.7) 

	}
}

genes<-list.files(dir,pattern='gene',full.names=T)
genes<-genes[grep('csv',genes)]
names<-gsub(basename(genes),pattern="_.*",replacement='')

for(i in 1:length(snps))
{
	phenotype<-info[info[,1,]%in%names[i],2]
	file<-read.table(genes[i],header=T,sep=',')
	dat<-subset(file,file$UncorrectedPvalue<=pval|file$CorrectedPvalue<=pval)[,1:3]
	if(nrow(dat)>0)
	{
	if(nrow(dat)>20)dat<-dat[1:20,]
	rownames(dat)<-1:nrow(dat)
	filt.xtable<-xtable(dat,caption=paste(phenotype,"Gene Results") ,digits=2, display = c(rep("s",2),rep("E",2)))
	print.xtable(filt.xtable, type="latex",file=paste0(dir,names[i],"_gene.tex"),scalebox=.7) 
	}
}



