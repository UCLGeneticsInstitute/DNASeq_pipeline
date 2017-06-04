message('Loading libraries')
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("GenomicRanges"))
suppressPackageStartupMessages(library("Gviz"))
#options(ucscChromosomeNames=FALSE)

option_list <- list(
	make_option("--skat", default=NULL,help="skat csv"),
	make_option("--sampleBam", default=NULL,help="optional bam file if reads are to be plotted",type='character'),
 	make_option("--outPDF",  help="pdf",type='character'), 
 	make_option("--Gene",  help="gene symbol eg ABCA4",type='character'),
 	make_option("--CNV",  help="CNV csv",type='character',default=NULL),
 	make_option("--Sample",  help="name of sample if track of samples variants to be included",type='character',default=NULL)
 	)

opt <- parse_args(OptionParser(option_list=option_list))
message("Starting Argument checks...") 
gen<-'hg19'

filt<-read.csv(opt$skat)
oPDF<-opt$outPDF

snp.file<-gsub(opt$skat,pattern='skat.csv',replacement='SKAT_results_by_SNP.tab')
if(!file.exists(snp.file)) stop('skat snp file doesnt exist')

test.bam<-opt$sampleBam
if(!is.null(test.bam))if(!file.exists(test.bam)) stop('bam doesnt exist')
CNV<-opt$CNV

gene<-opt$Gene
gene.row<-grep(gene,filt$Symbol)
if(length(gene.row)==0)stop('gene symbol not found in skat file')

sample<-opt$Sample
message('Finished argument check. processing.')

snp<-read.table(snp.file,header=T)
colnames(snp)<-c("SNP", "case.snp.hets", "case.snp.homs", "case.mafs.snp",
		"ctrl.snp.hets", "ctrl.snp.homs", "ctrl.mafs.snp", "ENSEMBL",
		"Symbol", "SKATO", "nb.snps", "nb.cases", "nb.ctrls", "nb.alleles.cases",
		"nb.alleles.ctrls", "case.maf", "ctrl.maf", "total.maf", "nb.case.homs",
		"nb.case.hets", "nb.ctrl.homs", "nb.ctrl.hets", "Chr", "Start",#
		"End", "FisherPvalue", "OddsRatio", "CompoundHetPvalue", "minCadd",
		"maxExac", "min.depth", "MeanCallRateCases", "MeanCallRateCtrls",
		"MaxMissRate", "HWEp", "MinSNPs", "MaxCtrlMAF", "SNPs", "GeneRD",
		"CaseSNPs", "SKATbeSNPs")
snp<-snp[snp$Symbol %in% filt$Symbol[gene.row],]
snp.data<-data.frame(t(data.frame(strsplit(snp$SNP,'_'))))
snp.data[,2]<-as.numeric(snp.data[,2])
rownames(snp.data)<-snp$SNP
colnames(snp.data)<-c('chr','start','a','b')
snp.data$case.maf<-as.numeric(snp$case.maf)
snp.data$case.maf[is.na(snp.data$case.maf)]<-0
snp.data$case.maf<- log10(as.numeric(snp.data$case.maf)+0.0001)
snp.data$ctrl.maf<-as.numeric(snp$ctrl.mafs.snp)
snp.data$ctrl.maf[is.na(snp.data$ctrl.maf)]<-0
snp.data$ctrl.maf<- log10(as.numeric(snp.data$ctrl.maf)+0.0001)

print(snp.data)
if(length(grep('chr',snp.data$chr))==0) snp.data$chr<-paste0('chr',snp.data$chr)
#snp.data$diff.maf<-snp.data$case.maf-snp.data$ctrl.maf
save(list=ls(environment()),file='prep.RData')

snp.data.ranges = GRanges(seqnames=snp.data$chr, ranges=IRanges(start=snp.data$start, end=snp.data$start) ,
	case=snp.data$case.maf,ctrl=snp.data$ctrl.maf)
colours<-c('red','black')
maf.data<- DataTrack(snp.data.ranges,col=colours,genome=gen,name='log10(MAF)',groups=c('cases','controls'))## case and control mafs

#snp.data.ranges = GRanges(seqnames=snp.data$chr, ranges=IRanges(start=snp.data$start, end=snp.data$start) ,
#	case=snp.data$diff.maf)
#maf.data<- DataTrack(snp.data.ranges,genome=gen,name='log10(MAF)',type='polygon')## case and control mafs

case.snps<-read.table(gsub(opt$skat,pattern='skat.csv',replacement='case_carriers'),header=FALSE)
case.snp.data<-case.snps[case.snps[,1] %in% rownames(snp.data),]
case.snp.data<-data.frame(cbind(case.snp.data,data.frame(t(data.frame(strsplit(case.snp.data[,1],'_'))))))
colnames(case.snp.data)<-c('SNP','Sample','ENSEMBL','HUGO','MinorAlleleCount','Chr','Start','Ref','Alt')
case.snp.data$Start<-as.numeric(case.snp.data$Start) ## all cases

if(!is.null(sample)) 
{
	sampleSNPs<-case.snp.data[case.snp.data$Sample %in% sample,]
	case.snp.ranges <- GRanges(seqnames=sampleSNPs$Chr, ranges=IRanges(start=sampleSNPs$Start, end=sampleSNPs$Start) ,
		case=sampleSNPs$MinorAlleleCount)
	target.snps<- DataTrack(case.snp.ranges,col='Red',genome=gen,name='CaseAlleleCount',ylim=range(sampleSNPs$MinorAlleleCount))
}

test.chr<-as.numeric(as.character(case.snp.data$Chr))
test.start<-min(as.numeric(as.character(case.snp.data$Start)) )
test.end<-max(as.numeric(as.character(case.snp.data$Start) )) 
test.dat<-unique(data.frame(chr=test.chr,start=test.start,end=test.end)) 
print(test.dat)

#dTrack4 <- DataTrack(range=test.bam, genome=gen, type="l", name="Coverage", window=-1, chromosome=test.chr)

itrack <- IdeogramTrack(genome = gen, chromosome = test.chr)

message('reading in exon bed file')
exon.bed<-read.table('/SAN/vyplab/UCLex/support/exons.hg19.bed')
bed.ranges = with(exon.bed, GRanges(V1, IRanges(start=V2, end=V3)))
test.dat$start<-as.numeric(test.dat$start)
test.dat$end<-as.numeric(test.dat$end)
test.ranges = with(test.dat, GRanges(chr, IRanges(start=start, end=end)))
hits<-data.frame( findOverlaps(bed.ranges ,test.ranges))
exon.bed.small<-exon.bed[hits$queryHits,]
colnames(exon.bed.small)<-c('chromosome','start','end','symbol')

exon.bed.small$width<-exon.bed.small$end-exon.bed.small$start
exon.bed.small$strand<-'+'
exon.bed.small$symbol<-gsub(exon.bed.small$symbol,pattern='_.*',replacement='') 
exon.bed.small$transcript<-exon.bed.small$symbol 
print(paste(nrow(exon.bed.small),'exons in gene')) 

genes<-unique(exon.bed.small$transcript)
exon.bed.small$exon<-NA
for(i in 1:length(genes))
{
	hit<-grep(genes[i],exon.bed.small$transcript)
	exon.bed.small$exon[hit]<-paste0(genes[i],'_',1:length(hit))
}

grtrack <- GeneRegionTrack(exon.bed.small, genome = gen, chromosome = test.chr, name = "Gene Model")
gtrack <- GenomeAxisTrack()

buffer<-1000


oData<-paste0(dirname(oPDF),'/plot.prep.RData')
save(list=ls(environment()),file=oData)
message(oData)


##### CNV
# add option to include exomeDepth call
if(!is.null(CNV))
{
	message('Plotting CNV call')
	# load CNV spreadsheet
	# extract sample and gene
	# plot track perhaps instead of bam file below. 
}

if(!is.null(test.bam)) 
{
	message('including bam track')
	bam.track <- AnnotationTrack(test.bam, genome = gen, chromosome = test.chr, name = "SampleBam")#ideally alignments track but it doesnt seem to work right now 
	pdf(oPDF,width=10,height=10)

	if(!is.null(sample)) 
	{
	plotTracks(list(itrack,gtrack, grtrack, bam.track,target.snps,maf.data), from=test.start, 
		to=test.end,chromosome=test.chr,extend.left=buffer,extend.right=buffer)
	}else #no sample specified so ignore that track 
	{
	plotTracks(list(itrack,gtrack, grtrack, bam.track,maf.data), from=test.start, 
		to=test.end,chromosome=test.chr,extend.left=buffer,extend.right=buffer)
	}
	dev.off()
	message(oPDF)
} else 
{
	pdf(oPDF)

	if(!is.null(sample)) 
	{
	plotTracks(list(itrack, gtrack, grtrack,target.snps,maf.data), from=test.start, 
		to=test.end,chromosome=test.chr,extend.left=buffer,extend.right=buffer)
	}else #no sample specified so ignore that track 
	{
	plotTracks(list(itrack, gtrack, grtrack,maf.data), from=test.start, 
		to=test.end,chromosome=test.chr,extend.left=buffer,extend.right=buffer)
	}
	dev.off()
	message(oPDF)
}
