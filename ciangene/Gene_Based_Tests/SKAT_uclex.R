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

######################
library(SKAT)

## these are teh snps that have been filterd by func and maf
snp.sp<-read.table(paste0(rootODir,'gene_test_variants_out.sp'),header=F)
snp.bim<-read.table(paste0(rootODir,'gene_test_variants_out.bim'),header=F)
snp.fam<-read.table(paste0(rootODir,'gene_test_variants_out.fam'),header=F)
colnames(snp.sp)<-snp.fam[,1]
rownames(snp.sp)<-snp.bim[,2]


## Three measures to choose valid controls I need to implement
#  ethnicity matching - remove non Caucasians
#  phenotype matching - removeConflictingControls in phenotype file making step - augment with matched control table
#	read depth matching - remove the samples who have poor coverage 


## ancestry matching
ancestry<-read.table(paste0(rootODir,'UCLex_samples_ancestry'),header=T)
caucasians<-ancestry$V1[ancestry$Caucasian]
snps<-snp.sp[,colnames(snp.sp)%in%caucasians]



## read depth 
read.depth<-read.table(paste0(rootODir,'Average.read.depth.by.gene.by.sample.tab'),header=T)
min.depth<-20
sample.means<-colMeans(read.depth[,4:ncol(read.depth)])
good.samples<-which(sample.means>=min.depth)

gene.means<-rowMeans(read.depth[,4:ncol(read.depth)])
good.genes<-which(gene.means>=min.depth)

snp.gene<-read.table(paste0(rootODir,'snp_gene_id'),header=F)

good.genes.snps<-snp.gene[snp.gene[,2] %in% read.depth[good.genes,1],1]


clean.snp.data<- snps[ rownames(snps) %in% good.genes.snps ,  colnames(snps) %in% colnames(read.depth)[good.samples] ]
## this is now a snp * sample matrix of snps that is filtered by
## gatk vqsr PASS 
## function and maf (refer to 'make_variant_list_gene_tests.R' for exact metrics )
## ancestry - non caucasians removed - may need to alter for IBD samples
## read depth - only the genes and samples that have a mean read depth by gene of 'min.depth' are kept. 



## Pheno matching
pheno.matching<-data.frame(read.csv('/SAN/vyplab/UCLex/scripts/DNASeq_pipeline/ciangene/Support/phenotype.info.csv',header=FALSE))


pheno<-read.table(paste0(rootODir,'Phenotypes'),header=F)
pheno.cohorts<-read.table(paste0(rootODir,'cohort.summary'),header=T)

good.genes.data<-snp.gene[snp.gene[,2] %in% read.depth[good.genes,1],]
uniq.genes<-unique(good.genes.data[,2])
nb.genes<-length(uniq.genes)


oDir<-paste0(rootODir,'SKAT/')
dir.create(oDir)



for(phen in 1:nrow(pheno.matching))
{
#	cases<-grep(pheno.matching[phen,1],pheno[,1])

	current.pheno<- grep(pheno.matching[phen,1],pheno.cohorts[,4])+2 ## plus two because first two rows of pheno file are ID

	remove.phenos<-length(which(pheno.matching[phen,3:ncol(pheno.matching)]>0))
	if(remove.phenos>0)
	{
		for(ct in 1:remove.phenos)
		{
			pheno[  grep(pheno.matching[phen,(ct+2)] , pheno[,1]),  current.pheno ] <-NA ## remove conflicting controls
		}
	}

	case.controls<-pheno[!is.na(pheno[,current.pheno]),1]
	clean.pheno.snps <- clean.snp.data[,colnames(clean.snp.data)%in%case.controls ]
	cases<-grep(pheno.matching[phen,1],colnames(clean.pheno.snps))

	cols<-c("Gene",'SKAT','SKATO','nb.snps')
	results<-data.frame(matrix(nrow=nb.genes,ncol=length(cols))) ## will add more columns later. (maf, call rate, mean depth etc)
	results[,1]<-uniq.genes
	results.out <- paste0(oDir,pheno.matching[phen,1],'_skat.csv')

	current.pheno<- rep(0,ncol(clean.pheno.snps))
	current.pheno[cases]<-1

	for(gene in 1:nb.genes)
	{
		gene.snps<-good.genes.data[ grep(uniq.genes[gene],good.genes.data[,2]) ,1]

		gene.snp.data<-clean.pheno.snps[ rownames(clean.pheno.snps) %in% gene.snps ,]
		nb.snps.in.gene<-nrow(gene.snp.data)
		print(paste(nb.snps.in.gene,'snps in', uniq.genes[gene]))

		if(nb.snps.in.gene)>0)
		{
			obj<-SKAT_Null_Model(current.pheno ~ 1, out_type="D")
			results[gene,2] <- SKAT(t(as.matrix(gene.snp.data)) , obj, missing_cutoff=0.4, estimate_MAF=2, kernel = "linear")$p.value
			results[gene,3] <- SKAT(t(as.matrix(gene.snp.data)) , obj, missing_cutoff=0.4, estimate_MAF=2, kernel = "linear")$p.value
			results[gene,3] <- nb.snps.in.gene
			write.table(data.frame(rownames(gene.snp.data),uniq.genes[gene]), paste0(rootODir,pheno.matching[phen,1],'_snps'), col.names=F,row.names=F,quote=F,sep='\t',append=T)
		}
	}

	results <- results[order(results[,2]), ]
	write.table(results,results.out,col.names=T,row.names=F,quote=F,sep='\t')
}