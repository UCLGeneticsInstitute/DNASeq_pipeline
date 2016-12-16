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
 ## changes: removed polyphen and sift fitlers.

 ## Firstly, I want ot make a list of variants to use for gene based tests

bim<-read.table(paste0(rootODir,'allChr_snpStats_out.bim'),header=F)
bim$SNP<-paste(bim[,1],bim[,4],bim[,5],sep="_")

dir<-paste0('/SAN/vyplab/UCLex/mainset_',release,'/mainset_',release, '_VEP/')
files <- list.files(dir, pattern ="csv$", full.names=T)
files <- files[order(as.numeric(gsub(gsub(basename(files), pattern ="chr", replacement =""), pattern = "_.*", replacement = "") ) )]

### filters
params<-read.table('/SAN/vyplab/UCLex/scripts/DNASeq_pipeline/ciangene/Support/params',header=FALSE)
min.maf<-params[grep('min.maf',params[,1]),2] # belwo this to be removed
max.maf<-params[grep('max.maf',params[,1]),2]  # above this to be removed
cadd<-params[grep('cadd',params[,1]),2] #below this to be removed

oFile<-paste0(rootODir,'filtered_snps_skat')
if(file.exists(oFile)) file.remove(oFile)

geneTable<-paste0(rootODir,'gene_dict_skat')
if(file.exists(geneTable)) file.remove(geneTable)

for(i in 1:length(files))
{
	system(paste("sed -i 's/#//g' ",files[i]))
	filt<-read.csv(files[i],header=T)
	filt$ExAC_MAF[is.na(filt$ExAC_MAF)]<-0
	filt.func.rare<-subset(filt,filt$ExAC_MAF<=max.maf & filt$ExAC_MAF>=min.maf & filt$CADD>= cadd)
	filt.func.rare$SNP<-paste(gsub(filt.func.rare$Location,pattern=':',replacement="_"),filt.func.rare$Allele,sep="_")
	filt.bim<-merge(filt.func.rare,bim,by='SNP')
	print(nrow(filt.bim))
	write.table(filt.bim$V2,oFile,col.names=FALSE,row.names=F,quote=F,sep='\t',append=TRUE)
	genes<-unique(data.frame(filt$Gene,filt$Feature,filt$SYMBOL))
	colnames(genes)<-c('ENSEMBL','ENST','Symbol')
	write.table(genes,geneTable,col.names=!file.exists(geneTable),row.names=F,quote=F,sep='\t',append=TRUE)
}