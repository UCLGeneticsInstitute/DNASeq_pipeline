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
library(S4Vectors) 
library(IRanges,lib.loc="/cluster/project8/vyp/cian/support/R") 
library(Rsamtools,lib.loc="/cluster/project8/vyp/cian/support/R")
library(Biostrings,lib.loc="/cluster/project8/vyp/cian/support/R")  
library(XVector,lib.loc="/cluster/project8/vyp/cian/support/R") 
library(GenomicRanges,lib.loc="/cluster/project8/vyp/cian/support/R") 
library(ExomeDepth,lib.loc="/cluster/project8/vyp/cian/support/R") 
data(exons.hg19)
data(Conrad.hg19)

fasta<-"/cluster/scratch3/vyp-scratch2/reference_datasets/human_reference_sequence/human_g1k_v37.fasta"
bamDir<-"/SAN/vyplab/UCLex/BAM/"

min.nb.per.batch<-5
max.nb.per.batch<-15

bams<-list.files(bamDir,pattern='bam$',full.names=TRUE)
bams<-bams[file.exists(bams)]
bamFiles<-data.frame(bams,basename(bams))
colnames(bamFiles)<-c('bampath','bam')
batches<-gsub(bamFiles$bam,pattern='_.*',replacement='')
batch.summary<-data.frame(table(batches))
batch.summary<-batch.summary[order(-batch.summary$Freq),]

bamFiles$batch<-batches
colnames(bamFiles)<-c('bampath','bam','batch')
for(batch in 1:nrow(batch.summary))
{
	batch.samples<-grep(batch.summary$batches[batch],bamFiles$batch)
	if(batch.summary$Freq[batch]>max.nb.per.batch) ## if the gruops are too big i want to split them up
	{
		group.batches<-split(bamFiles$bam[batch.samples], ceiling(seq_along(bamFiles$bam[batch.samples])/max.nb.per.batch))
		for(lis in 1:length(group.batches))
		{
			subgroup<-group.batches[[lis]]
			subgroup.hit<- bamFiles$bam %in% subgroup
			bamFiles$batch[subgroup.hit]<-paste0(bamFiles$batch[subgroup.hit],'_',lis)
		}
	}
	if(batch.summary$Freq[batch]<min.nb.per.batch)bamFiles$batch[batch.samples]<-'Remove'
}



outDir<-paste0(rootODir,'ExomeDepth/')
if(!file.exists(outDir))dir.create(outDir)

write.table(bamFiles,paste0(outDir,'Sample.Batch.Guide.tab'),col.names=T,row.names=F,quote=F,sep='\t')
bamFiles<-bamFiles[!bamFiles$batch=='Remove',]
uniq.batches<-unique(bamFiles$batch)


for(batch in 1:length(uniq.batches))
{
	batch_outDir<-paste0(outDir,uniq.batches[batch],'/')
	if(!file.exists(batch_outDir))dir.create(batch_outDir)
	countsFile<-paste0(batch_outDir,"counts.RData")

	batch.bams<-bamFiles$bampath[bamFiles$batch %in% uniq.batches[batch]]
	
	if(length(batch.bams)>1)
	{
		if(!file.exists(countsFile))
		{
			print("Getting counts per region") 
			#my.counts <- getBamCounts(bed.file='bins_for_exomeDepth.bed',
			my.counts <- getBamCounts(bed.frame = exons.hg19,,
				bam.files = batch.bams,
				include.chr = FALSE,
				referenceFasta = fasta)
			# turn counts into dataframe
			my.counts.dat <- as(my.counts[, colnames(my.counts)], 'data.frame')
			print(head(my.counts.dat))
			my.counts.dat$chromosome <- gsub(as.character(my.counts.dat$space),pattern = 'chr',replacement = '')
			print(uniq.batches[batch])
			save.image(file=countsFile ) 
		} else load(countsFile)
	}
}

exit

for(case in 1:length(samples))
{
	my.test <- my.counts.dat[,grep(basename(samples[case]),colnames(my.counts.dat) )  ]
	my.ref.samples <- samples[-case]
	my.reference.set <- as.matrix(my.counts.dat[, basename(my.ref.samples) ])
	my.choice <- select.reference.set (test.counts = my.test,
		reference.counts = my.reference.set,
		bin.length = (my.counts.dat$end - my.counts.dat$start)/1000,
		n.bins.reduced = 10000)
	my.matrix <- as.matrix( my.counts.dat[, my.choice$reference.choice, drop = FALSE])
	my.reference.selected <- apply(X = my.matrix,
		MAR = 1,
		FUN = sum)
	all.exons <- new('ExomeDepth',
		test  =my.test,
		reference = my.reference.selected,
		formula = 'cbind(test, reference) ~ 1')

	all.exons <- CallCNVs(x = all.exons,
		transition.probability = 10^-4,
		chromosome = my.counts.dat$space,
		start = my.counts.dat$start,
		end = my.counts.dat$end,
		name = my.counts.dat$names)
	if(nrow(all.exons@CNV.calls)>0)
	{
		all.exons <- AnnotateExtra(x = all.exons,
		reference.annotation = Conrad.hg19.common.CNVs,
		min.overlap = 0.5,
		column.name = 'Conrad.hg19')

		exons.hg19.GRanges <- GenomicRanges::GRanges(seqnames = exons.hg19$chromosome,
		IRanges::IRanges(start=exons.hg19$start,end=exons.hg19$end),
		names = exons.hg19$name)
		
		all.exons <- AnnotateExtra(x = all.exons,
		reference.annotation = exons.hg19.GRanges,
		min.overlap = 0.0001,
		column.name = 'exons.hg19')

		callFile<-paste0(outDir,gsub(basename(samples[case]),pattern='\\.bam',replacement=''),'_calls.RData') 
		save(all.exons,file=callFile)
		output.file <- paste0(outDir,gsub(basename(samples[case]),pattern='\\.bam',replacement=''),'_calls.csv') 
		write.csv(file = output.file, x = all.exons@CNV.calls,row.names = FALSE)
	} 
}
