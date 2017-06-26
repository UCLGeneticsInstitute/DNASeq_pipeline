library(ExomeDepth)
bayes.filter<-0

data.directory<-'/SAN/vyplab/UCLex/mainset_July2016/cian/ExomeDepth/'
data.files<-list.files(data.directory,pattern='calls.csv',recursive=TRUE,full.names=TRUE)
data.names<-gsub(basename(data.files),pattern='csv',replacement='RData')

for(file in 1:length(data.files))
{
	dat<-read.csv(data.files[file])
	dat$Sample<-data.names[file]
	if(file==1)allCalls<-dat else allCalls<-data.frame(rbind(allCalls,dat)) 
}

genes<-read.table('/SAN/vyplab/UCLex/support/genes.bed',header=TRUE,sep='\t')

allCalls$Gene<-gsub(allCalls$exons.hg19,pattern='_.*',replacement='')
allCalls<-allCalls[order(allCalls$Gene),]
allCalls$batch<-substr(allCalls$Sample,1,11)
allCalls$batch<-gsub(allCalls$batch,pattern='b',replacement='B') 

allCalls$short.name<-gsub(allCalls$Sample,pattern='_sor.*',replacement='')
pdf.out<-paste0(data.directory,'CNVplots.pdf')
pdf(pdf.out)

mac.genes<-read.table('/cluster/project8/vyp/cian/data/Support/CandidateGenes/candidate_genes_macular_dystrophy')[,1]
allCalls<-allCalls[allCalls$Gene %in% mac.genes,]

out.bed<-paste0(data.directory,'TargetGenes.bed') 
if(file.exists(out.bed)) file.remove(out.bed)

for(cnv in 1:nrow(allCalls))
{
	cnv.file<-paste0(data.directory,allCalls$batch[cnv],'/',allCalls$Sample[cnv])
	if(file.exists(cnv.file))
	{
		print(cnv.file)
		load(cnv.file)
		flank<-1000
		calls<- all.exons@CNV.calls
		calls<-subset(calls,calls$BF>=bayes.filter)
		calls$Gene<-gsub(calls$exons.hg19,pattern='_.*',replacement='')

		current.cnv<-calls[calls$Gene %in% allCalls$Gene[cnv],]
		cnv.id<-paste0(allCalls$short.name[cnv][1],'_',current.cnv$Gene[1])

		gene.start<-min(current.cnv$start)-flank
		gene.end<-max(current.cnv$end)+flank
		chr<-unique(current.cnv$chromosome) 

		if(cnv==1)cnvs.plotted<-'Dud'
		if(length(grep(cnv.id,cnvs.plotted))==0) # im plotting once per CNV containing gene per sample so skip CNVs after first one. 
		{
			plot (all.exons,
			sequence = chr,
			xlim = c(gene.start, gene.end),
			count.threshold = 20,
			main = paste(allCalls$short.name[cnv][1],current.cnv$Gene[1]) ,
			cex.lab = 0.8,
			with.gene = FALSE)

			gene.bed<-data.frame(chr,gene.start,gene.end) ## record which gnees ahve CNVs in them. will use BEd for conifer afterwards.
			gene.bed$id<-paste(unlist(gene.bed),collapse='_')
			if(length(gene.bed$id)>1) stop('IDs for multiple CNVs not merged properly')
			if(cnv==1)all.bed<-gene.bed
			if(cnv>1)
			{
				if(length( grep(gene.bed$id,all.bed$id) )== 0 ) all.bed<-data.frame(rbind(all.bed,gene.bed))
			} 
		}
		if(cnv==1)cnvs.plotted<-cnv.id else cnvs.plotted<-data.frame(rbind(cnvs.plotted,cnv.id))

	}
}
dev.off()

write.table(all.bed,out.bed,col.names=F,row.names=F,quote=F,sep='\t',append=FALSE)
message(paste('BED file made in',out.bed))
