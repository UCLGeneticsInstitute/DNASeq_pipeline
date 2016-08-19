getArgs <- function() {
  myargs.list <- strsplit(grep("=",gsub("--","",commandArgs()),value=TRUE),"=")
  myargs <- lapply(myargs.list,function(x) x[2] )
  names(myargs) <- lapply(myargs.list,function(x) x[1])
  return (myargs)
}

rootODir<-'/SAN/vyplab/UCLex/mainset_July2016/cian/'
release<-'July2016'

myArgs <- getArgs()

if ('rootODir' %in% names(myArgs))  rootODir <- myArgs[[ "rootODir" ]]
if ('release' %in% names(myArgs))  release <- myArgs[[ "release" ]]

#####################
#library(SKAT)
#library(parallel)
#library(snpStats)


## these are teh snps that have been filterd by func and maf
snp.sp<-read.table(paste0(rootODir,'gene_test_variants_out.sp'),header=F)
snp.bim<-read.table(paste0(rootODir,'gene_test_variants_out.bim'),header=F)
snp.fam<-read.table(paste0(rootODir,'gene_test_variants_out.fam'),header=F)
colnames(snp.sp)<-snp.fam[,1]
rownames(snp.sp)<-snp.bim[,2]


snp.gene.base<-read.table(paste0(rootODir,'snp_gene_id'),header=F)
colnames(snp.gene.base)<-c("SNP",'ENSEMBL')
gene.dict<-read.table(paste0(rootODir,'gene_dict_skat'),header=T)
gene.dict$ENST<-NULL
gene.dict<-unique(gene.dict)

snp.gene<-merge(gene.dict,snp.gene.base,by="ENSEMBL",all.y=T)

oDir<-paste0(rootODir,'SKATnew/')
if(!file.exists(oDir))dir.create(oDir)

qc<-paste0(oDir,'/qc/')
if(!file.exists(qc))dir.create(qc)

test<-FALSE


## Three measures to choose valid controls I need to implement
#  ethnicity matching - remove non Caucasians
#  phenotype matching - removeConflictingControls in phenotype file making step - augment with matched control table
#	read depth matching - remove the samples who have poor coverage 


## ancestry matching
ancestry<-read.table(paste0(rootODir,'UCLex_samples_ancestry'),header=T)
caucasians<-ancestry$V1[ancestry$Caucasian]
snps<-snp.sp[,colnames(snp.sp)%in%caucasians]



## read depth 
min.depth<-20

read.depth<-read.table(paste0(rootODir,'Average.read.depth.by.gene.by.sample.tab'),header=T)

sample.means<-colMeans(read.depth[,5:ncol(read.depth)])
good.samples<-which(sample.means>=min.depth)
bad.samples<-names(which(sample.means<min.depth))
write.table(bad.samples,paste0(qc,'samples_removed_because_of_low_read_depth.tab'),col.names=F,row.names=F,quote=F,sep='\t')

gene.means<-rowMeans(read.depth[,5:ncol(read.depth)])
good.genes<-which(gene.means>=min.depth)
bad.genes<-read.depth$Gene[which(gene.means<min.depth)]
write.table(bad.genes,paste0(qc,'genes_removed_because_of_low_read_depth.tab'),col.names=F,row.names=F,quote=F,sep='\t')


good.genes.snps<-snp.gene$SNP[snp.gene$ENSEMBL %in% read.depth$Gene[good.genes]]
clean.snp.data<- snps[ rownames(snps) %in% good.genes.snps ,  colnames(snps) %in% colnames(read.depth)[good.samples] ]
snp.gene.clean<-snp.gene[snp.gene$SNP %in%rownames(clean.snp.data),]
write.table(snp.gene.clean,paste0(qc,'read.depth.ancestry.func.clean.snps.tab'),col.names=F,row.names=F,quote=F,sep='\t')

## this is now a snp * sample matrix of snps that is filtered by
## gatk vqsr PASS 
## function and maf (refer to 'make_variant_list_gene_tests.R' for exact metrics )
## ancestry - non caucasians removed - may need to alter for IBD samples
## read depth - only the genes and samples that have a mean read depth by gene of 'min.depth' are kept. 



## Pheno matching
pheno.matching<-data.frame(read.csv('/SAN/vyplab/UCLex/scripts/DNASeq_pipeline/ciangene/Support/phenotype.info.csv',header=FALSE))


pheno<-read.table(paste0(rootODir,'Phenotypes'),header=F)
pheno.cohorts<-read.table(paste0(rootODir,'cohort.summary'),header=T)

good.genes.data<-snp.gene[snp.gene$ENSEMBL %in% read.depth$Gene[good.genes],]
uniq.genes<-unique(good.genes.data$ENSEMBL)
nb.genes<-length(uniq.genes)
exit

##for(phen in 1:nrow(pheno.matching)) skip the first few samples. 
for(phen in 70:nrow(pheno.matching))
{
#	cases<-grep(pheno.matching[phen,1],pheno[,1])

	current.pheno<- grep(paste0(pheno.matching[phen,1],'$'),pheno.cohorts[,4])+2 ## plus two because first two rows of pheno file are ID
	## $ added at end of grep to differentiate between batches Lambiase and LambiaseSD. 

	remove.phenos<-length(which(pheno.matching[phen,3:ncol(pheno.matching)]>0))
	if(remove.phenos>0)
	{
		for(ct in 1:remove.phenos)
		{
			pheno[  grep(pheno.matching[phen,(ct+2)] , pheno[,1]),  current.pheno ] <-NA ## remove conflicting controls
		}
	}

	ex.ctrl.in<-paste0(rootODir,'External_Control_data/', pheno.matching[phen,1], '_ex_ctrls.frq_clean')
	if(file.exists(ex.ctrl.in))## When preceding scripts are fully run, this file should always exist.  
	{
	ex.ctrl<-read.table(ex.ctrl.in,header=T)
	ex.ctrl.rare.snps<-ex.ctrl$SNP[ex.ctrl$MAF <= 0.01]
	} else ex.ctrl.rare.snps<-rownames(clean.snp.data) 

	case.controls<-pheno[!is.na(pheno[,current.pheno]),1]
	clean.pheno.snps <- clean.snp.data[rownames(clean.snp.data) %in% ex.ctrl.rare.snps ,colnames(clean.snp.data)%in%case.controls ]
	write.table(rownames(clean.snp.data), paste0(oDir,pheno.matching[phen,1],'_filt_snp_list'), col.names=F,row.names=F,quote=F,sep='\t')
	cases<-grep(pheno.matching[phen,1],colnames(clean.pheno.snps))

	if(length(cases)>20)
	{
	cols<-c("Gene",'SKATO','nb.snps','nb.cases','nb.ctrls')
	results<-data.frame(matrix(nrow=nb.genes,ncol=length(cols))) ## will add more columns later. (maf, call rate, mean depth etc)
	colnames(results)<-cols
	results$Gene<-uniq.genes

	results<-merge(gene.dict,results,by.y='Gene',by.x='ENSEMBL',all.y=T)
	srt<-data.frame(1:length(uniq.genes),uniq.genes)
	results<-merge(results,srt,by.y='uniq.genes',by.x='ENSEMBL')
	results<-results[order(results[,(ncol(results))]),]
	results<-results[,1:(ncol(results)-1)]

	uniq.genes<-unique(results[,1])
	nb.genes<-length(uniq.genes)

	current.pheno<- rep(0,ncol(clean.pheno.snps))
	current.pheno[cases]<-1
	results$nb.cases<-table(current.pheno)[[2]]
	results$nb.ctrls<-table(current.pheno)[[1]]

	## to speed up, run 4 tests per time. theres probably a better way to do this. 
	oData<-paste0(oDir,pheno.matching[phen,1],'.RData')
	save(snp.gene,nb.genes,good.genes.data,uniq.genes,clean.pheno.snps,current.pheno,results,phen,pheno.matching,oDir , file=oData )
	script.out<-paste0(oDir,pheno.matching[phen,1],'.R')
	oData<-paste0('load("',paste0(oDir,pheno.matching[phen,1],'.RData"'),')')
	write.table(oData,script.out,col.names=F,row.names=F,quote=F,sep='\t')
	runR='sh /cluster/project8/vyp/cian/scripts/bash/runRonCluster.sh'
	file.append( script.out, 'SKAT_uclex.template.R') 
	system(paste(runR,script.out))


	if(test)
	{
		 ## this bit moved to SKAT_uclex.template.R for now - makes it easier to run more phenotypes at once on cluster. 
	}# test 
		
}

} 
exit
files<-list.files(oDir,pattern="sh",full.names=T)
mclapply(files,function(x)system(paste("sh",x)),mc.cores=4)

