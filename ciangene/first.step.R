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
library(snpStats)

params<-read.table('/SAN/vyplab/UCLex/scripts/DNASeq_pipeline/ciangene/Support/params',header=FALSE)
percent.ext.ctrls <-params[grep('min.maf',params[,1]),2] /100

dir<-paste0('/SAN/vyplab/UCLex/mainset_',release,'/mainset_',release, '_snpStats/') 

files <- list.files(dir, pattern ="_snpStats.RData", full.names=T) 
files <- files[order(as.numeric(gsub(gsub(basename(files), pattern ="chr", replacement =""), pattern = "_.*", replacement = "") ) )]
print(files)

############################################## New - combined annovar with VEP 
files.chr<-as.numeric( gsub(gsub(basename(files),pattern='chr',replacement=''),pattern='_.*',replacement=''))

func <- c("nonsynonymous SNV", "stopgain SNV", "nonframeshift insertion", "nonframeshift deletion", "frameshift deletion",
		"frameshift substitution", "frameshift insertion",  "nonframeshift substitution", "stoploss SNV", "splicing"
		,"exonic;splicing")
lof <-  c("frameshift deletion", "frameshift substitution", "frameshift insertion",  "stoploss SNV", "splicing"
		,"stopgain SNV","exonic;splicing"
		)

dir<-paste0('/SAN/vyplab/UCLex/mainset_',release,'/mainset_',release, '_VEP/') 
vep.list <- list.files(dir, pattern ="csv$", full.names=T) 
vep.chr<-as.numeric( gsub(gsub(basename(vep.list),pattern='.*chr',replacement=''),pattern='\\..*',replacement=''))

params<-read.table('/SAN/vyplab/UCLex/scripts/DNASeq_pipeline/ciangene/Support/params',header=FALSE)
min.maf<-params[grep('min.maf',params[,1]),2] # belwo this to be removed
max.maf<-params[grep('max.maf',params[,1]),2]  # above this to be removed
cadd<-params[grep('cadd',params[,1]),2] #below this to be removed


#oFile<-paste0(rootODir,'Annotations/annovar.vep.merged.anno.tab')
#if(file.exists(oFile))file.remove(oFile)
func.file<-paste0(rootODir,'Annotations/func.tab')
if(file.exists(func.file))file.remove(func.file)
lof.file<-paste0(rootODir,'Annotations/lof.tab')
if(file.exists(lof.file))file.remove(lof.file)
############################################




print(files)

oDir <- rootODir

readDepthDir <- file.path(rootODir,"read_depth")

if(!file.exists(readDepthDir) ) dir.create(readDepthDir)
readDepthoFile <- "Depth_Matrix.sp"

full <- file.path(oDir, "allChr_snpStats.sp") ## added '.sp' suffix
annotations.dir<-file.path(oDir, "Annotations") 
if(!file.exists(annotations.dir) ) dir.create(annotations.dir)
annotations.out <- paste0(annotations.dir,"annotations.snpStat")

oMap <- paste0(oDir, "/UCLex_", release, ".map")
oBim <- paste0(oDir, "/UCLex_", release, ".bim")

FIRST <- TRUE

for(i in 1:length(files)){
#for (i in c(22)) {
	message("Now loading file ", files[i])
	load(files[i])

##      Extract clean variants. 
	pass.snps <- annotations.snpStats$clean.signature[which(annotations.snpStats$FILTER == "PASS")] #
	matrix.calls.snpStats <- matrix.calls.snpStats[, which( colnames(matrix.calls.snpStats) %in% pass.snps ) ]#
	annotations.snpStats <- subset(annotations.snpStats, annotations.snpStats$FILTER == "PASS")
	matrix.depth <-  matrix.depth[which( rownames(matrix.depth) %in% pass.snps) ,]
	matrix.depth <- apply(matrix.depth, 2, as.numeric)

	# Make map file 
	map <- data.frame(matrix(nrow=nrow(annotations.snpStats), ncol = 4) ) 
	map[,1]<-annotations.snpStats$Chr
	map[,2]<-rownames(annotations.snpStats) 
	map[,3] <- 0 
	map[,4]<-annotations.snpStats$Start

	matrix.calls.snpStats<-matrix.calls.snpStats[,order(map[,4])]
	annotations.snpStats<-annotations.snpStats[order(map[,4]),]
	matrix.depth<-matrix.depth[order(map[,4]),]
	map<-map[order(map[,4]),]
	if(FIRST) 
	{	
		samples <- rownames(matrix.calls.snpStats)
		if(percent.ext.ctrls>0)
		{
			ext.ctrls <- sample(length(samples), length(samples) * percent.ext.ctrls) # Get list of extCtrls. 
			write.table(length(ext.ctrls) , file = paste0(oDir, "nb_ext_ctrl_samples"), col.names=F, row.names=F, quote=F, sep="\t") 
			write.table(rownames(matrix.calls.snpStats[ext.ctrls,]) , file = paste0(oDir, "ext_ctrl_samples"), col.names=F, row.names=F, quote=F, sep="\t") 
	  		ext.samples <- matrix.calls.snpStats[ext.ctrls ,]
			ext.samples.sum <- data.frame(colnames(matrix.calls.snpStats), col.summary(ext.samples) ) 
			ext.samples.names <- data.frame(rownames(ext.samples) , row.summary(ext.samples) ) 
			write.table(ext.samples.sum, file = paste0(oDir, "Ext_ctrl_variant_summary") , col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t", append=FALSE) 
			write.table(ext.samples.names, file = paste0(oDir, "Ext_ctrl_sample_summary") , col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t", append=FALSE)
		}

		write.SnpMatrix(t(matrix.calls.snpStats), full, col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t", append=FALSE)   ##this is where it becomes numbers
   		write.table(annotations.snpStats, annotations.out, col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t", append=FALSE)
		write.table(matrix.depth, file = file.path(readDepthDir, readDepthoFile) , col.names=FALSE, row.names=FALSE, quote = FALSE, sep="\t" , append = FALSE)
		write.table(data.frame('SNP'=colnames(matrix.calls.snpStats), col.summary(matrix.calls.snpStats)),paste0(oDir,"UCLex_",release,"_all_samples_variant_summary"),col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")

		write.table(map, oMap, col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t", append=FALSE) 
		write.table(data.frame(cbind(map,annotations.snpStats$Obs,annotations.snpStats$Ref ) ) , oBim, col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t", append=FALSE)
        FIRST <- FALSE
  	} else {
  		if(percent.ext.ctrls>0)
		{
	  		ext.samples <- matrix.calls.snpStats[ext.ctrls ,]
			ext.samples.sum <- data.frame(colnames(matrix.calls.snpStats), col.summary(ext.samples) ) 
			ext.samples.names <- data.frame(rownames(ext.samples) , row.summary(ext.samples) ) 
	   		write.table(ext.samples.sum, file = paste0(oDir, "Ext_ctrl_variant_summary") , col.names=F, row.names=F, quote=F, sep="\t", append=T) 
			write.table(ext.samples.names, file = paste0(oDir, "Ext_ctrl_sample_summary") , col.names=F, row.names=F, quote=F, sep="\t", append=T) 
		}
		write.SnpMatrix(t(matrix.calls.snpStats), full, col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t", append=TRUE)  
   		write.table(annotations.snpStats, annotations.out, col.names=FALSE, row.names=TRUE, quote=FALSE, sep="\t", append=TRUE) 
		write.table(matrix.depth,  file = file.path(readDepthDir, readDepthoFile) , col.names=FALSE, row.names=FALSE, quote = FALSE, sep="\t" , append = TRUE) 
		write.table(data.frame('SNP'=colnames(matrix.calls.snpStats), col.summary(matrix.calls.snpStats)),paste0(oDir,"UCLex_",release,"_all_samples_variant_summary"),col.names=F,row.names=F,quote=F,sep="\t",append=T) 

		write.table(map, oMap, col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t", append=TRUE)  
		write.table(data.frame(cbind(map, annotations.snpStats$Obs,annotations.snpStats$Ref ) ) , oBim, col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t", append=TRUE) 
	}

############################################ New
	vep<-read.csv(vep.list[vep.chr %in% files.chr[i]],header=T)
	vep$ExAC_MAF[is.na(vep$ExAC_MAF)]<-0
	snplist<-data.frame(rownames(annotations.snpStats),annotations.snpStats$ExonicFunc)
	colnames(snplist)<-c('SNP','ExonicFunc')
	id<-data.frame(t(data.frame(strsplit(snplist$SNP,'_'))))
	snplist$id<-paste0(id[,1],':',id[,2])
	merged.snp.data<-merge(snplist,vep,by.x='id',by.y='Location')

	merged.snp.data.keep <-data.frame(merged.snp.data$SNP,merged.snp.data$ExonicFunc,merged.snp.data$Existing_variation,merged.snp.data$Gene,merged.snp.data$Consequence,merged.snp.data$SYMBOL,merged.snp.data$ExAC_MAF,merged.snp.data$CADD)
	colnames(merged.snp.data.keep)<-gsub(colnames(merged.snp.data.keep),pattern='.*\\.',replacement='')

	merged.snp.data.keep$Func<-FALSE
	merged.snp.data.keep$LOF<-FALSE
	merged.snp.data.keep$Rare<-FALSE
	merged.snp.data.keep$Rare[merged.snp.data.keep$ExAC_MAF <= max.maf] <- TRUE
	merged.snp.data.keep$Func[merged.snp.data.keep$ExonicFunc %in% func] <- TRUE
	merged.snp.data.keep$LOF[merged.snp.data.keep$ExonicFunc %in% lof] <- TRUE
	merged.rare<-subset(merged.snp.data.keep,merged.snp.data.keep$Rare)

	merged.func<-subset(merged.rare,merged.rare$Func)
	write.table(merged.func,func.file,col.names=!file.exists(func.file),row.names=FALSE,quote=FALSE,sep='\t',append=TRUE)
	merged.func<-subset(merged.rare,merged.rare$LOF)
	write.table(merged.func,lof.file,col.names=!file.exists(lof.file),row.names=FALSE,quote=FALSE,sep='\t',append=TRUE)
############################################

} ## for(i in 1:length(files)

## Make fam file. 
fam <- data.frame(matrix(nrow=ncol(matrix.depth), ncol = 6)  ) 
fam[,1] <- colnames(matrix.depth)
fam[,2] <- colnames(matrix.depth)
fam[,3] <- 0 
fam[,4] <- 0 
fam[,5] <- 1 
fam[,6] <- 0
write.table(fam, paste0(oDir, "/UCLex_", release, ".fam") , col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t") 


#bim <- read.table(oBim, header=FALSE, sep="\t") 
#bim <- bim[with(bim, order(V1, V4)), ]
#write.table(bim, oBim, col.names=FALSE, row.names=FALSE, quote=FALSE ,sep="\t") 
