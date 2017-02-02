loopPlot<-function(dat,CallsDirectory)
{
	for(cnv in 1:nrow(dat))
	{
		cnv.file<-paste0(CallsDirectory,'/single_exons/',dat$sample[cnv],'.RData')
		if(!grepl('CNVcalls',cnv.file))cnv.file<-paste0(CallsDirectory,'/CNVcalls/single_exons/',dat$sample[cnv],'.RData')
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
				affected.genes<-unique(genes[hits$queryHits,]$Name)
				aff.bed<-genes[genes$Name%in% affected.genes,]
				gene.start<-min(aff.bed$start)-flank
				gene.end<-max(aff.bed$end)+flank

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
