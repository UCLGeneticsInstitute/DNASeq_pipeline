suppressPackageStartupMessages(library(ExomeDepth) )
suppressPackageStartupMessages(library(optparse) )
suppressPackageStartupMessages(library(GenomicRanges) )
suppressPackageStartupMessages(library(biomaRt) )

dbup='/cluster/project8/vyp/cian/scripts/bash/dropbox_uploader.sh'
source('/cluster/project8/vyp/cian/scripts/r/ExomeDepthplot.R')
option_list <- list(
	make_option(c("-v", "--verbose"), action="store_true", default=TRUE,help="Print extra output [default]"),
 	make_option(c("--CallsDirectory"),  help="location of calls from exomedepth",type='character',default=NULL),
 	make_option(c("--outPDF"),  help="where to store pdf",type='character',default=NULL),
 	make_option(c("--Genes"),  help="candidate gene list",type='character',default=NULL),
 	make_option(c("--BF"),  help="candidate gene list",type='character',default=NULL),
 	make_option(c("--novelBF"),  help="candidate gene list",type='character',default=NULL),
 	make_option(c("--Pattern"),  default='bam.cnv',type='character'),
 	make_option(c("--SavePrep"), default=TRUE, help="Do you want to save an image of setup?",type='character'),
 	make_option(c("--TargetChr"),  default=NULL,type='character',help='Are there certain chromosomes of interest?'),
 	make_option(c("--DBOXdir"),  default=NULL,type='character')
 )


opt <- parse_args(OptionParser(option_list=option_list))
if ( opt$verbose ) {
 write("Starting Argument checks...\n", stderr())
}
######################
CallsDirectory<-opt$CallsDirectory
if(!file.exists(CallsDirectory))dir.create(CallsDirectory)
outPDF<-opt$outPDF
candidate.genes<-opt$Genes
pattern<-opt$Pattern
bayes.filter<-as.numeric(opt$BF)
novel.bayes.filter<-opt$novelBF
DBOXdir<-opt$DBOXdir
SavePrep<-opt$SavePrep
if(!is.null(opt$TargetChr))TargetChrs<-paste(unlist(strsplit(opt$TargetChr,',') ))  else TargetChrs<-opt$TargetChr

message(paste('Getting calls from:',CallsDirectory))

data.directory<-paste0(dirname(outPDF),'/') 
data.files<-list.files(CallsDirectory,pattern=pattern,recursive=TRUE,full.names=TRUE)
data.names<-gsub(data.files,pattern=pattern,replacement='bam_X.RData')
data.names<-gsub(data.names,pattern='multi_exons',replacement='single_exons')

for(file in 1:length(data.files))
{
	dat<-read.csv(data.files[file],sep='\t')
	if(file==1)allCalls<-dat else allCalls<-data.frame(rbind(allCalls,dat)) 
}

ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", host="grch37.ensembl.org")
genes <- getBM(attributes=c('ensembl_gene_id','gene_biotype','hgnc_symbol','chromosome_name','start_position','end_position'), filters = 'chromosome_name', values ="X", mart = ensembl)
gr0 = with(genes, GRanges(chromosome_name, IRanges(start=start_position, end=end_position)))

allCalls<-allCalls[order(allCalls$genePos_hg19.tab),]
allCalls$short.name<-gsub(allCalls$sample,pattern='_sor.*',replacement='')

print(paste('Number of CNVs found across all genes:',nrow(allCalls)))

#######################################################################################################################################
## function to plot all CNVs listed in data.frame
loopPlot<-function(dat)
{
	for(cnv in 1:nrow(dat))
	{
		cnv.file<-paste0(CallsDirectory,'single_exons/',dat$sample[cnv],'_X.RData')
		if(file.exists(cnv.file))
		{
			load(cnv.file)
			flank<-1000
			calls<- sample.mod@CNV.calls
			current.cnv<-calls[calls$id %in% dat$id[cnv],]
			cnv.id<-paste0(dat$sample[cnv],'_',current.cnv$id[1])
			print(cnv.id)

			chr<-unique(current.cnv$chromosome) 
			gr1 <- with(current.cnv, GRanges(chromosome, IRanges(start=start, end=end)))
			hits<-data.frame(findOverlaps(gr0, gr1)) 
			
			if(nrow(hits)>0) ## If CNV is in gene(s) in bed file then use same Xaxis for each CNV to make plots prettier
			{
				affected.genes<-unique(genes[hits$queryHits,]$hgnc_symbol)
				aff.bed<-genes[genes$hgnc_symbol%in% affected.genes,]
				gene.start<-min(aff.bed$start_position)-flank
				gene.end<-max(aff.bed$end_position)+flank

			} else
			{
				gene.start<-min(current.cnv$start)-flank
				gene.end<-max(current.cnv$end)+flank
			}

			if(cnv==1)cnvs.plotted<-'Dud'
			if(length(grep(cnv.id,cnvs.plotted))==0) # im plotting once per CNV containing gene per sample so skip CNVs after first one. 
			{
				message(paste("CNV file",cnv.file))
				pl(sample.mod,
				sequence = chr,
				xlim = c(gene.start, gene.end),
				count.threshold = 20,
				main = paste(dat$sample[cnv],dat$id[cnv]), 
				cex.lab = 0.8,
				with.gene = TRUE)

				gene.bed<-data.frame(chr,gene.start,gene.end) ## record which gnees ahve CNVs in them. will use BEd for conifer afterwards.
				gene.bed$id<-paste(unlist(gene.bed),collapse='_')
				if(length(gene.bed$id)>1) stop('IDs for multiple CNVs not merged properly')

			} else message(cnv.id)
			if(cnv==1)cnvs.plotted<-cnv.id else cnvs.plotted<-data.frame(rbind(cnvs.plotted,cnv.id))
		} else message(paste(cnv.file, 'doesnt exist...'))
	}
}

if(SavePrep)
{
	robj<-paste0(data.directory,'test_setup.RData')
	message(paste('Saving workspace image to', robj))
	save(list=ls(environment()),file=robj)
}	
#######################################################################################################################################
#######################################################################################################################################
##  Plot the most significnat novel multi exon CNVs. 
novelCNVs<-read.table(list.files(CallsDirectory,pattern='multi_exons_postQC_CNVs_X.tab',recursive=TRUE,full.names=TRUE),header=TRUE)
novelCNVs<-novelCNVs[with(novelCNVs, order(chromosome,start,end)), ]


allCNVsPDF<-paste0(dirname(outPDF),'/AllCNVs.pdf')
pdf(allCNVsPDF)
	loopPlot(novelCNVs) 
dev.off()

gene.CNVs<- data.frame(table(unlist(strsplit(novelCNVs$genePos_hg19.tab,',' ) ) ))
gene.CNVs<-gene.CNVs[order(-gene.CNVs$Freq),]
write.table(gene.CNVs,paste0(dirname(outPDF),'/GeneCNVcount.csv'),col.names=T,row.names=F,quote=T,sep=',')

if(!is.null(novel.bayes.filter))novelCNVs.sig<-subset(novelCNVs,novelCNVs$BF >= novel.bayes.filter)else novelCNVs.sig<-novelCNVs
write.table(novelCNVs.sig,paste0(dirname(outPDF),'/NovelCNVs.csv'),col.names=T,row.names=F,quote=T,sep=',')
message(paste(nrow(novelCNVs.sig), 'novel post-QC CNVs found that pass Bayes Filter'))


plot.singletons<-FALSE
if(plot.singletons)
{
	cnv.freq<-data.frame(table(novelCNVs.sig$id)) 
	novelCNVs.sig$unique<-FALSE
	for(i in 1:nrow(novelCNVs.sig))
	{
		hit.row<-grep(as.character(novelCNVs.sig$id[i]) ,as.character(cnv.freq[,1]))
		if(cnv.freq[hit.row,2]==1) novelCNVs.sig$unique[i]<-TRUE
	}
	print( table(novelCNVs.sig$unique))
	novelCNVs.sig<-subset(novelCNVs.sig,novelCNVs.sig$unique)

	#out.bed<-paste0(data.directory,'NovelGenes.bed') 
	#if(file.exists(out.bed)) file.remove(out.bed)
	novelCNVspdf<-paste0(dirname(outPDF),'/NovelCNVs.pdf')

	pdf(novelCNVspdf)
	loopPlot(novelCNVs.sig)
	dev.off()
	message(paste('PDF is:',novelCNVspdf))
}



hitCNVspdf<-paste0(dirname(outPDF),'/hitCNVs.pdf')
cnv.genes<-subset(gene.CNVs,gene.CNVs$Freq>1) # Plot CNVs in genes with more than one CNV
novelCNVs$interesting<-FALSE
for(cn in 1:length(cnv.genes)){
	hit.rows<-grep(cnv.genes[cn,1],novelCNVs$genePos_hg19.tab)
	novelCNVs$interesting[hit.rows]<-TRUE
}
topCNVs<-subset(novelCNVs,novelCNVs$interesting & novelCNVs$Conrad_CNVs_hg19.tab=='none')

topCNVs$CNVid<-paste0(topCNVs$chromosome,'-',topCNVs$sample)
topCNVs$CNVplot<-1:nrow(topCNVs)
max.distance<-10000
for(cnv in 1:nrow(topCNVs)) # i want to merge rows that are basically the same CNV
{
	cnvs<-grep(topCNVs$CNVid[cnv],topCNVs$CNVid) # so here I find the CNVs that are nearby
	if(length(cnvs)>1)
	{
		for(cn in 2:length(cnvs))
		{
			if( topCNVs[cnvs,]$start[cn] - topCNVs[cnvs,]$end[cn-1] < max.distance ) 
			{
				topCNVs$CNVplot[topCNVs$CNVid == topCNVs[cnvs,]$CNVid[cn] ] <- topCNVs[cnvs,]$CNVplot[cn]
			}
		}
	}
}
cnvs<-unique(topCNVs$CNVplot)
for(cnv in 1:length(cnvs)) # and here I merge them. 
{
	topCNVs$start[topCNVs$CNVplot==cnvs[cnv]]<- min(topCNVs$start[topCNVs$CNVplot==cnvs[cnv]])
	topCNVs$end[topCNVs$CNVplot==cnvs[cnv]]<- max(topCNVs$end[topCNVs$CNVplot==cnvs[cnv]])
}
topCNVs.merged<-topCNVs[unique(topCNVs$CNVplot),]
topCNVs.merged<-topCNVs.merged[with(topCNVs.merged, order(chromosome, start,end)), ]
pdf(hitCNVspdf)
loopPlot(topCNVs.merged)
dev.off()
#if(file.info(hitCNVspdf) $size<10000) file.remove(hitCNVspdf) else message(paste('PDF is:',hitCNVspdf))
#######################################################################################################################################
## Now look at the candidate genes. 

if(!is.null(candidate.genes))
{
	mac.genes<-read.table(candidate.genes)[,1]
	allCalls$CandidateGene<-FALSE
	allCalls$CandidateGeneName<-NA

	for(i in 1:length(mac.genes))
	{
		hit<-grep(mac.genes[i],allCalls$genePos_hg19.tab,ignore.case=TRUE)
		if(length(hit)>0)print(hit)
		if(length(hit)>0)allCalls$CandidateGene[hit]<-TRUE
		if(length(hit)>0)allCalls$CandidateGeneName[hit]<-mac.genes[i]
	}
	allCalls.candidate<-allCalls[allCalls$CandidateGene,]

	if(nrow(allCalls.candidate)==0)
	{
		message("No CNVs in candidate genes :(")
	} else 
	{
		data.directory<-paste0(dirname(outPDF),'/') 
		#out.bed<-paste0(data.directory,'TargetGenes.bed') 
		#if(file.exists(out.bed)) file.remove(out.bed)
		if(!is.null(bayes.filter))allCalls.candidate<-subset(allCalls.candidate,allCalls.candidate$BF >= bayes.filter)

		pdf(outPDF)
		loopPlot(allCalls.candidate)
		dev.off()
		if(file.info(outPDF) $size<10000) file.remove(outPDF) else message(paste('PDF is:',outPDF))

	}
}

################################################################################################
### Plot target chromosomes
############################################################################################
if(!is.null(TargetChrs))
{
	message(paste('Plotting Target Chromosomes:',TargetChrs))
	chr.cnvs<-novelCNVs[novelCNVs$chromosome %in% TargetChrs,]
	targetChrPDF<-paste0(dirname(outPDF),'/TargetChrsCNVs.pdf')

	message(paste(nrow(chr.cnvs),'CNVs on target chromosomes... plotting in',targetChrPDF)) 
	pdf(targetChrPDF)
		loopPlot(chr.cnvs)
	dev.off()
}


################################################################################################
### Push plots to Dropbox
############################################################################################
if(!is.null(DBOXdir))
{
	DBOXdir<-paste0(DBOXdir,'/')
	run<-paste(dbup,'upload',outPDF, paste0('PostDoc/', DBOXdir, basename(outPDF))) 
	system(run)
	run<-paste(dbup,'upload',allCNVsPDF, paste0('PostDoc/', DBOXdir, basename(allCNVsPDF))) 
	system(run)
	if(plot.singletons)
	{	
		run<-paste(dbup,'upload',novelCNVs, paste0('PostDoc/', DBOXdir,basename(novelCNVspdf))) 
		system(run)
	}
	run<-paste(dbup,'upload',paste0(dirname(outPDF),'/NovelCNVs.csv'), paste0('PostDoc/', DBOXdir,'NovelCNVs.csv')) 
	system(run)
	run<-paste(dbup,'upload',paste0(dirname(outPDF),'/GeneCNVcount.csv'), paste0('PostDoc/', DBOXdir,'GeneCNVcount.csv')) 
	system(run)
	run<-paste(dbup,'upload',hitCNVspdf, paste0('PostDoc/', DBOXdir,'genesMultipleCNVs.pdf')) 
	if(file.exists(hitCNVspdf))system(run)
	if(!is.null(TargetChrs))run<-paste(dbup,'upload',targetChrPDF, paste0('PostDoc/', DBOXdir,basename(targetChrPDF) ) ) 
	if(!is.null(TargetChrs))system(run)

}


message("Finished ok")