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

dir<-paste0('/SAN/vyplab/UCLex/mainset_',release,'/mainset_',release, '_snpStats/') 

files <- list.files(dir, pattern ="_snpStats.RData", full.names=T) 
files <- files[order(as.numeric(gsub(gsub(basename(files), pattern ="chr", replacement =""), pattern = "_.*", replacement = "") ) )]

print(files)

sample.cohort<-read.table(paste0(rootODir,'Sample.cohort'),header=F)
uniq.cohorts<-unique(sample.cohort[,2])
nb.cohorts<-length(uniq.cohorts)

FIRST<-TRUE
for(i in 1:length(files))
{
	load(files[i])
	print(files[i])
	matrix.depth <- apply(matrix.depth, 2, as.numeric)


	for(long in 1:length(annotations.snpStats$Gene)) ## Some variant names label with multiple genes, so i want to separate and include all of tehse
	{
		gene<-annotations.snpStats$Gene[long]
		if(nchar(gene)>15)
		{
		#	long.hit<- gregexpr(pattern ='ENS', x)
		#	for(i in 1:length(long.hit[[1]]))
		#	{
		#		genes<-c(genes, substr(x,long.hit[[1]][i],long.hit[[1]][i]+14 ) )
		#	}
			annotations.snpStats$Gene[long]<-substr(gene,1,15)
		} else annotations.snpStats$Gene[long]<-gene
	}

	uniq.genes<-unique(annotations.snpStats$Gene)
	uniq.genes.dat<-data.frame(matrix(nrow=length(uniq.genes),ncol=4+nb.cohorts,0))
	colnames(uniq.genes.dat)<-c('Gene','Chr','Start','End',uniq.cohorts)
	uniq.genes.dat$Gene<-uniq.genes

	uniq.genes.sample<-data.frame(matrix(nrow=length(uniq.genes),ncol=4+ncol(matrix.depth),0))
	colnames(uniq.genes.sample)<-c('Gene','Chr','Start','End',colnames(matrix.depth))
	uniq.genes.sample$Gene<-uniq.genes




	for(gene in 1:length(uniq.genes))
	{
		hits<-which(annotations.snpStats$Gene==uniq.genes[gene])
		gene.start.row<-min(hits)
		gene.end.row<-max(hits)
		chr<-unique(annotations.snpStats$Chr[gene.start.row:gene.end.row])
		start<-min(annotations.snpStats$Start[gene.start.row:gene.end.row])
		end<-max(annotations.snpStats$End[gene.start.row:gene.end.row])
		
		uniq.genes.dat$Chr[gene]<-chr
		uniq.genes.dat$Start[gene]<-start
		uniq.genes.dat$End[gene]<-end
		for(batch in 1:nb.cohorts)
		{
			cohort.hits<-sample.cohort[,2]%in%uniq.cohorts[batch]
			mean.depth<-mean(matrix.depth[gene.start.row:gene.end.row,cohort.hits],na.rm=T) # average readDepth per gene by cohort
			uniq.genes.dat[gene,4+batch]<-mean.depth
		}


		uniq.genes.sample$Chr[gene]<-chr
		uniq.genes.sample$Start[gene]<-start
		uniq.genes.sample$End[gene]<-end
		for(sample in 1:ncol(matrix.depth))
		{
			mean.depth<-mean(matrix.depth[gene.start.row:gene.end.row,sample],na.rm=T) # average readDepth per gene by sample
			uniq.genes.sample[gene,4+batch]<-mean.depth
		}
	}#gene
	print(nrow(uniq.genes.dat))
	if(FIRST)
	{
		write.table(uniq.genes.dat,paste0(rootODir,'Average.read.depth.by.gene.by.cohort.tab'),col.names=T,row.names=F,quote=F,sep='\t')
		write.table(uniq.genes.sample,paste0(rootODir,'Average.read.depth.by.gene.by.sample.tab'),col.names=T,row.names=F,quote=F,sep='\t')
		FIRST<-FALSE
	} else
	{
		write.table(uniq.genes.dat,paste0(rootODir,'Average.read.depth.by.gene.by.cohort.tab'),col.names=F,row.names=F,quote=F,sep='\t',append=T)
		write.table(uniq.genes.sample,paste0(rootODir,'Average.read.depth.by.gene.by.sample.tab'),col.names=F,row.names=F,quote=F,sep='\t',append=T)
	}
}