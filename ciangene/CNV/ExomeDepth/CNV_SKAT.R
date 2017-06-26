suppressPackageStartupMessages(library(optparse) )
dbup='/cluster/project8/vyp/cian/scripts/bash/dropbox_uploader.sh'
option_list <- list(
	make_option(c("--CNVdir"),  type='character',help="Print extra output [default]"),
 	make_option(c("--SKATdir"),  help="skat output dir",type='character',default=NULL),
 	make_option(c("--outFile"),  help="oDir",type='character',default=NULL)
)
opt <- parse_args(OptionParser(option_list=option_list))
write("Starting Argument checks...\n", stderr())

######################
CNVdir<-opt$CNVdir
SKATdir<-opt$SKATdir
outFile<-opt$outFile
CNVfiles<-list.files(CNVdir,recursive=TRUE,pattern='NovelCNVs.csv',full.names=TRUE)
print(CNVfiles)
gene.cnv.count<-list.files(CNVdir,recursive=TRUE,pattern='GeneCNVcount.csv',full.names=TRUE)
if(!file.exists(dirname(outFile)))dir.create(CallsDirectory,recursive=TRUE)

for(i in 1:length(gene.cnv.count))
{
	dat<-data.frame(batch=basename(dirname(gene.cnv.count[i])),read.csv(gene.cnv.count[i])) 
	if(i==1)gene.count<-dat else gene.count<-data.frame(rbind(gene.count,dat))
}
gene.flt<-gene.count[-grep('male',gene.count$batch),]
gene.means<-aggregate(gene.flt$Freq,by=list(gene.flt$Var1),mean)
cnv.laden.genes<-subset(gene.means,gene.means$x>5)

for(i in 1:length(CNVfiles))
{
	dat<-data.frame(batch=basename(dirname(CNVfiles[i])),read.csv(CNVfiles[i])) 
	if(length(grep('male',dat$batch))>0)
	{
		dat<-subset(dat,dat$chromosome>22)
		if(nrow(dat)>0)dat<-data.frame(dat, "multi_exons_postQC_CNVs.tab"  ,   "ID")
	}
	if(nrow(dat)>0)
	{
		if(length(grep('IDs',colnames(dat)))>0) dat$IDs<-NULL
		if(i!=1)colnames(dat)<-colnames(cnvs)
		if(length(grep('cnvs',ls()))==0)cnvs<-dat else cnvs<-data.frame(rbind(cnvs,dat))
	}
}


skat.files<-list.files(SKATdir,recursive=TRUE,full.names=TRUE)
cases<-read.csv(skat.files[grep('case.list',skat.files )],header=FALSE)
solved<-read.csv(skat.files[grep('solved',skat.files )],header=FALSE)
cases$SKATgenes<-NA
cases$SKATvariants<-NA
cases$CNVgenes<-NA
cases$CNVid<-NA
for(i in 1:nrow(cases))
{
	case.skat<-solved[grep(cases[i,1],solved$Case),]
	if(nrow(case.skat)>1)
	{
		cases$SKATgenes[i]<- paste(unique(unlist(strsplit(case.skat$MostLikelyCausativeGenesInOrder,';') )),collapse=';')
		cases$SKATvariants[i]<- paste(unique(unlist(strsplit(case.skat$Variants,';') )),collapse=';')
	}
	case.cnv<-cnvs[grep(cases[i,1],cnvs$sample),]
	if(nrow(case.cnv)>0)
	{
		cases$CNVgenes[i]<- paste(unlist(case.cnv$ID),collapse=';')
		cases$CNVid[i]<- paste(unlist(case.cnv$id),collapse=';')
	}
}

write.table(cases,outFile,col.names=T,row.names=F,quote=T,sep=',')
out<-paste('summary file written to', outFile)
message(out)