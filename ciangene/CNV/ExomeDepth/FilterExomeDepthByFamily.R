suppressPackageStartupMessages(library(ExomeDepth) )
suppressPackageStartupMessages(library(optparse) )
suppressPackageStartupMessages(library(data.table) )
suppressPackageStartupMessages(library(GenomicRanges) )


dbup='/cluster/project8/vyp/cian/scripts/bash/dropbox_uploader.sh'
source('/cluster/project8/vyp/cian/scripts/r/ExomeDepthplot.R')
option_list <- list(
 	make_option("--inputFile",  help='iFile',default=NULL,type='character'),
 	make_option("--proband",  help="pband",type='character',default=NULL),
 	make_option("--family",  help="fam",type='character',default=NULL),
 	make_option("--outDir",  help="odir",type='character',default=NULL)
 )

opt <- parse_args(OptionParser(option_list=option_list))
message('starting argument check')
inputFile<-opt$inputFile
message(paste('input file:',inputFile)) 
proband<-opt$proband
family<-opt$family
outDir<-opt$outDir
if(!file.exists(outDir))dir.create(outDir,recursive=TRUE)

message('Loading plot functions')
source('/SAN/vyplab/UCLex/scripts/DNASeq_pipeline/ciangene/CNV/ExomeDepth/loopPlot.R')


file<-fread(inputFile)
proband<-fread(proband)
family<-fread(family)

proband.cnvs<-file[file$sample %in% proband,]
family.cnvs<-file[file$sample %in% family,]

proband.ranges <- with(proband.cnvs, GRanges(chromosome, IRanges(start=start, end=end)))
family.ranges <- with(family.cnvs, GRanges(chromosome, IRanges(start=start, end=end)))

hits<-data.frame(findOverlaps(proband.ranges, family.ranges)) 

proband.novel.cnvs<-proband.cnvs[-unique(hits$queryHits),]
if(nrow(proband.novel.cnvs)==0)
{
	stop('There are no novel CNVs in proband apparently')
} else
{
	write.table(proband.novel.cnvs,paste0(outDir,'ProbandOnlyCNVs.csv'),col.names=T,row.names=F,quote=T,sep=',')
	pdfOut<-paste0(outDir,'ProbandOnlyCNVs.pdf')
	pdf(pdfOut)
		loopPlot(proband.novel.cnvs)
	dev.off()
	message(paste('output in csv & pdf at',paste0(outDir,'ProbandOnlyCNVs.csv and .pdf') ))
}

