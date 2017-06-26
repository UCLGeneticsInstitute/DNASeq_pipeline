suppressPackageStartupMessages(library(S4Vectors) ) 
suppressPackageStartupMessages(library(IRanges) )
suppressPackageStartupMessages(library(Rsamtools))
suppressPackageStartupMessages(library(Biostrings))  
suppressPackageStartupMessages(library(XVector) )
suppressPackageStartupMessages(library(GenomicRanges)) 
suppressPackageStartupMessages(library(ExomeDepth) )
suppressPackageStartupMessages(library(optparse) )

option_list <- list(
	make_option(c("--chrom"), default=NULL,help="Chromosome"),
	make_option(c("-v", "--verbose"), action="store_true", default=TRUE,help="Print extra output [default]"),
 	make_option(c("--SampleList"),  help="one column file containing list of samples",type='character',default=NULL),
 	make_option(c("--BamList"),  help="one column file containing list of samples",type='character',default=NULL),
 	make_option(c("--oDir"), help="oDir",type='character')

 )


opt <- parse_args(OptionParser(option_list=option_list))
if ( opt$verbose ) {
 write("Starting Argument checks...\n", stderr())
}
######################
SampleList<-opt$SampleList
BamList<-opt$BamList


data(exons.hg19)
data(Conrad.hg19)

fasta<-"/cluster/scratch3/vyp-scratch2/reference_datasets/human_reference_sequence/human_g1k_v37.fasta"
min.nb.per.batch<-5
max.nb.per.batch<-15

if(!is.null(SampleList) & !is.null(BamList) ) stop("Specify either list of samples to be found or BAM file list directly. Not both.")

if(!is.null(SampleList))
{
	if(!file.exists(SampleList))stop("SampleList file doesn't exist. ")
	message(paste('Reading samples from',paste0('--',SampleList,'--')))
	SampleList<-read.table(SampleList,header=FALSE)
	message(paste('Provided',nrow(SampleList),'samples'))
	SampleList<-SampleList[,1]

	## Find bams for the samples
	bamList<-list.files('/SAN/vyplab/UCLex_raw',pattern='bam$',full.names=T) 
	bamList<-bamList[grep('sorted_unique',bamList)]
	sample.bams<-vector()
	sample.names<-vector()
	for(i in 1:length(SampleList))
	{
		hit<-grep(SampleList[i],bamList)
		if(length(hit)>0)sample.bams<-c(sample.bams,bamList[hit])
		if(length(hit)>0)sample.names<-c(sample.names,SampleList[i])
		if(length(hit)==0)print(paste('No BAM file found for',SampleList[i]))
	}
	message(paste('Found BAM files for',length(sample.bams),'samples'))

	if(length(SampleList)>max.nb.per.batch) message('There are too many samples, consider splitting?')
	if(length(SampleList)<min.nb.per.batch) message('There are too few samples, results may be unreliable')
	if(length(sample.bams)>length(SampleList))
	{
		message("Multiple BAMs exist for >=1 samples. Will use the larger one")
		#Find the duplicate files and remove the smaller one. 
		dups<-sample.bams[grep(sample.names[which( lapply(lapply(sample.names,function(x) grep(x,sample.bams)),function(x) length(x)) >1) ],sample.bams)]
		file.sizes<-file.size(dups)
		sample.bams<-sample.bams[sample.bams!=dups[file.sizes==min(file.sizes)]]
	}
} 

if(!is.null(BamList))
{
	sample.bams<-read.table(BamList,header=FALSE)[,1]
	message(paste('Read file containing',length(sample.bams),'sample BAMs'))
}

outDir<-paste0(opt$oDir,'/') 
if(!file.exists(outDir))dir.create(outDir)
message(paste('Results will be placed in --',outDir))

write.table(sample.bams,paste0(outDir,'SampleList'),col.names=F,row.names=F,quote=F,sep='\t')

countsFile<-paste0(outDir,"/counts.RData")
if(!file.exists(countsFile))
{
	print("Getting counts per region") 
	my.counts <- getBamCounts(bed.frame = exons.hg19,,
			bam.files = sample.bams,
			include.chr = FALSE,
			referenceFasta = fasta)
		# turn counts into dataframe
		my.counts.dat <- as(my.counts[, colnames(my.counts)], 'data.frame')
		print(head(my.counts.dat))
		my.counts.dat$chromosome <- gsub(as.character(my.counts.dat$space),pattern = 'chr',replacement = '')
		save.image(file=countsFile ) 
} else load(countsFile)

message('Finished getting counts. now going to call CNVs')

sample.cols<-grep('bam',colnames(my.counts.dat))
genes<-read.table('/SAN/vyplab/UCLex/support/genes.bed',header=TRUE,sep='\t')

for(case in 1:length(sample.cols))
{
	print(paste('Current case is:',colnames(my.counts.dat)[sample.cols][case]))
	output.file<-paste0(outDir,colnames(my.counts.dat)[sample.cols[case]], '_CNV_calls.csv')
	callFile<-paste0(outDir,colnames(my.counts.dat)[sample.cols[case]], '_CNV_calls.RData')
	pdfFile<-paste0(outDir,colnames(my.counts.dat)[sample.cols[case]], '_CNV_calls.pdf')
	#pdf(pdfFile)

	my.test <- my.counts.dat[,sample.cols[case] ]
	my.reference.set <- as.matrix(my.counts.dat[,sample.cols[-case] ])

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

		save(all.exons,file=callFile)
		write.csv(file = output.file, x = all.exons@CNV.calls,row.names = FALSE)
	} 
}