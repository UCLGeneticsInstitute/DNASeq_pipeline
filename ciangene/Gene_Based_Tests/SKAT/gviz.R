## plot sample snps as hom/het

message('Loading libraries')
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("GenomicRanges"))
suppressPackageStartupMessages(library("Gviz"))

option_list <- list(
	make_option("--skat", default=NULL,help="Chromosome"),
	make_option("--sampleBam", default=NULL,help="Print extra output [default]",type='character'),
 	make_option("--outPDF",  help="pdf",type='character'), 
 	make_option("--Gene",  help="pdf",type='character'),
 	make_option("--Sample",  help="pdf",type='character')
 	)

opt <- parse_args(OptionParser(option_list=option_list))
message("Starting Argument checks...") 
gen<-'hg19'

filt<-read.csv(opt$skat)
oPDF<-opt$outPDF

snp.file<-gsub(opt$skat,pattern='SKAT_processed_filtered.csv',replacement='results_by_SNP_cases.csv')
if(!file.exists(snp.file)) stop('skat snp file doesnt exist')

test.bam<-opt$sampleBam
if(!file.exists(test.bam)) stop('bam doesnt exist')

gene<-opt$Gene
gene.row<-grep(gene,filt$Symbol)
if(length(gene.row)==0)stop('gene symbol not found in skat file')

message('Finished argument check. processing.')

snp<-read.csv(snp.file)
snp<-snp[snp$Symbol %in% filt$Symbol[gene.row],]
snps<-unlist(strsplit(snp[,1],';'))
snp.data<-data.frame(t(data.frame(strsplit(snps,'_'))))
snp.data[,2]<-as.numeric(snp.data[,2])
rownames(snp.data)<- snps

colnames(snp.data)<-c('chr','start','a','b')
snp.data$case.maf<-snp$maf.snp.cases
snp.data$case.maf[is.na(snp.data$case.maf)]<-0
snp.data$case.maf<-log10(as.numeric(snp.data$case.maf)+0.0001)
snp.data$ctrl.maf<-snp$maf.snp.ctrls
snp.data$ctrl.maf[is.na(snp.data$ctrl.maf)]<-0
snp.data$ctrl.maf<-log10(as.numeric(snp.data$ctrl.maf)+0.0001)

case.snps<-gsub(opt$skat,pattern='SKAT_processed_filtered.csv',replacement='case_carriers.csv',header=FALSE)
case.snp.data<-case.snps[case.snps[,1] %in% rownames(snp.data),]
case.snp.data<-data.frame(cbind(case.snp.data,data.frame(t(data.frame(strsplit(case.snp.data[,1],'_'))))))
case.snp.ranges = GRanges(seqnames=snp.data$chr, ranges=IRanges(start=snp.data$start, end=snp.data$start) ,
	case=snp.data$case.maf)
target.snps<- DataTrack(snp.data.ranges,groups=colours,col=colours,genome=gen,name='log10(MAF)')


snp.data.ranges = GRanges(seqnames=snp.data$chr, ranges=IRanges(start=snp.data$start, end=snp.data$start) ,
	case=snp.data$case.maf,ctrl=snp.data$ctrl.maf)

test.chr<-filt$Chr[gene.row]
test.start<-filt$Start[gene.row]
test.end<-filt$End[gene.row]
test.dat<-data.frame(chr=test.chr,start=test.start,end=test.end)
print(test.dat)

dTrack4 <- DataTrack(range=test.bam, genome=gen, type="l", name="Coverage", window=-1, chromosome=test.chr)
itrack <- IdeogramTrack(genome = gen, chromosome = test.chr)

message('reading in exon bed file')
exon.bed<-read.table('/SAN/vyplab/UCLex/support/exons.hg19.bed')
bed.ranges = with(exon.bed, GRanges(V1, IRanges(start=V2, end=V3)))
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
#bmt <- BiomartGeneRegionTrack(genome="hg19", chromosome=test.chr, start=test.start, end=test.end,
                        #   filter=list(with_ox_refseq_mrna=TRUE), stacking="dense")
gtrack <- GenomeAxisTrack()


colours<-c('red','black')
maf.data<- DataTrack(snp.data.ranges,groups=colours,col=colours,genome=gen,name='log10(MAF)')





if(!is.null(test.bam)) 
{
	message('including bam track')
	alTrack <- AlignmentsTrack(test.bam, isPaired=TRUE,start=test.start,end=test.end,chromosome=test.chr,genome='hg19')
	pdf(oPDF)
	plotTracks(list(itrack, gtrack, alTrack,  grtrack,maf.data), from=test.start, 
		to=test.end,chromosome=test.chr,extend.left=1000,extend.right=1000)
	dev.off()
	message(oPDF)
} else 
{
	pdf(oPDF)
	plotTracks(list(itrack, gtrack,  grtrack,maf.data), from=test.start, 
		to=test.end,chromosome=test.chr,extend.left=1000,extend.right=1000)
	dev.off()
	message(oPDF)
}

save(list=ls(environment()),file='prep.RData')

