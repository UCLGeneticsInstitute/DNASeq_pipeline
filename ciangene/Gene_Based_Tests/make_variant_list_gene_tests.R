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
library(MultiPhen)
library(parallel)


## Firstly, I want ot make a list of variants to use for gene based tests

bim<-read.table(paste0(rootODir,'allChr_snpStats_out.bim'),header=F)
bim$SNP<-paste(bim[,1],bim[,4],bim[,5],sep="_")

dir<-paste0('/SAN/vyplab/UCLex/mainset_',release,'/mainset_',release, '_VEP/') 

files <- list.files(dir, pattern ="csv$", full.names=T) 
files <- files[order(as.numeric(gsub(gsub(basename(files), pattern ="chr", replacement =""), pattern = "_.*", replacement = "") ) )]



### filters
min.maf<-0 # belwo this to be removed
max.maf<-0.01 # above this to be removed
cadd<-15 #below this to be removed

polyphen.remove<-'benign'
sift.remove<-'tolerated'


#func <- c("nonsynonymous SNV", "stopgain SNV", "nonframeshift insertion", "nonframeshift deletion", "frameshift deletion",
	#	"frameshift substitution", "frameshift insertion",  "nonframeshift substitution", "stoploss SNV", "splicing"
	#	,"exonic;splicing")
#lof <-  c("frameshift deletion", "frameshift substitution", "frameshift insertion",  "stoploss SNV", "splicing"
	#	,"stopgain SNV","exonic;splicing"
	#	)

oFile<-paste0(rootODir,'filtered_snps_skat')

for(i in 1:length(files))
{
	system(paste("sed -i 's/#//g' ",files[i]))
	file<-read.csv(files[i],header=T)
	filt<-file[-grep(polyphen.remove,file$PolyPhen), ]
	filt<-filt[-grep(sift.remove,filt$SIFT), ]
	filt.func.rare<-subset(filt,filt$ExAC_MAF<=max.maf & filt$ExAC_MAF>=min.maf & filt$CADD>= cadd)
	filt.func.rare$SNP<-paste(gsub(filt.func.rare$Location,pattern=':',replacement="_"),filt.func.rare$Allele,sep="_")

	filt.bim<-merge(filt.func.rare,bim,by='SNP')
	print(nrow(filt.bim))
	if(i==1)
	{
		write.table(filt.bim$V2,oFile,col.names=F,row.names=F,quote=F,sep='\t')
	} else
	{
		write.table(filt.bim$V2,oFile,col.names=F,row.names=F,quote=F,sep='\t',append=T)
	}
}