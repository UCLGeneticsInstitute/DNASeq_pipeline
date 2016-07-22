getArgs <- function() {
  myargs.list <- strsplit(grep("=",gsub("--","",commandArgs()),value=TRUE),"=")
  myargs <- lapply(myargs.list,function(x) x[2] )
  names(myargs) <- lapply(myargs.list,function(x) x[1])
  return (myargs)
}

release <- 'September2015'
rootODir<-'/scratch2/vyp-scratch2/cian'

myArgs <- getArgs()

if ('rootODir' %in% names(myArgs))  rootODir <- myArgs[[ "rootODir" ]]
if ('release' %in% names(myArgs))  release <- myArgs[[ "release" ]]

######################
library(snpStats)
percent.ext.ctrls <- .10


dir <- paste0("/cluster/project8/vyp/exome_sequencing_multisamples/mainset/GATK/mainset_", release , "/mainset_", release,  "_snpStats/")
files <- list.files(dir, pattern ="_snpStats.RData", full.names=T) 
files <- files[order(as.numeric(gsub(gsub(basename(files), pattern ="chr", replacement =""), pattern = "_.*", replacement = "") ) )]

print(files)

oDir <- paste0(rootODir, "/UCLex_", release, "/")
if(!file.exists(oDir)) dir.create(oDir)

readDepthDir <- paste0(rootODir, "/UCLex_", release, "/read_depth/")
if(!file.exists(readDepthDir) ) dir.create(readDepthDir)
readDepthoFile <- "Depth_Matrix.sp"

full <- paste0(oDir, "/allChr_snpStats.sp") ## added '.sp' suffix
annotations.out <- paste0(oDir, "/annotations.snpStat")
out <- paste0(full, ".RData") 
a.out <- paste0(annotations.out, ".RData") 

oMap <- paste0(oDir, "/UCLex_", release, ".map")
oBim <- paste0(oDir, "/UCLex_", release, ".bim")


for(i in 1:length(files)){
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
	if(i==1) 
	{	
		samples <- rownames(matrix.calls.snpStats)
		ext.ctrls <- sample(length(samples), length(samples) * percent.ext.ctrls) # Get list of extCtrls. 
		write.table(length(ext.ctrls) , file = paste0(oDir, "nb_ext_ctrl_samples"), col.names=F, row.names=F, quote=F, sep="\t") 
		write.table(rownames(matrix.calls.snpStats[ext.ctrls,]) , file = paste0(oDir, "ext_ctrl_samples"), col.names=F, row.names=F, quote=F, sep="\t") 
  		ext.samples <- matrix.calls.snpStats[ext.ctrls ,]
		ext.samples.sum <- data.frame(colnames(matrix.calls.snpStats), col.summary(ext.samples) ) 
		ext.samples.names <- data.frame(rownames(ext.samples) , row.summary(ext.samples) ) 

		write.table(ext.samples.sum, file = paste0(oDir, "Ext_ctrl_variant_summary") , col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t", append=FALSE) 
		write.table(ext.samples.names, file = paste0(oDir, "Ext_ctrl_sample_summary") , col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t", append=FALSE)
		write.SnpMatrix(t(matrix.calls.snpStats), full, col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t", append=FALSE)   ##this is where it becomes numbers
   		write.table(annotations.snpStats, annotations.out, col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t", append=FALSE)
		write.table(matrix.depth, file = paste0(readDepthDir, "/", readDepthoFile) , col.names=FALSE, row.names=FALSE, quote = FALSE, sep="\t" , append = FALSE)
		write.table(data.frame('SNP'=colnames(matrix.calls.snpStats), col.summary(matrix.calls.snpStats)),paste0(oDir,"UCLex_",release,"_all_samples_variant_summary"),col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")

		write.table(map, oMap, col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t", append=FALSE) 
		write.table(data.frame(cbind(map,annotations.snpStats$Obs,annotations.snpStats$Ref ) ) , oBim, col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t", append=FALSE)
  	} else {
  		ext.samples <- matrix.calls.snpStats[ext.ctrls ,]
		ext.samples.sum <- data.frame(colnames(matrix.calls.snpStats), col.summary(ext.samples) ) 
		ext.samples.names <- data.frame(rownames(ext.samples) , row.summary(ext.samples) ) 
   		write.table(ext.samples.sum, file = paste0(oDir, "Ext_ctrl_variant_summary") , col.names=F, row.names=F, quote=F, sep="\t", append=T) 
		write.table(ext.samples.names, file = paste0(oDir, "Ext_ctrl_sample_summary") , col.names=F, row.names=F, quote=F, sep="\t", append=T) 

		write.SnpMatrix(t(matrix.calls.snpStats), full, col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t", append=TRUE)  
   		write.table(annotations.snpStats, annotations.out, col.names=FALSE, row.names=TRUE, quote=FALSE, sep="\t", append=TRUE) 
		write.table(matrix.depth,  file = paste0(readDepthDir, "/", readDepthoFile) , col.names=FALSE, row.names=FALSE, quote = FALSE, sep="\t" , append = TRUE) 
		write.table(data.frame('SNP'=colnames(matrix.calls.snpStats), col.summary(matrix.calls.snpStats)),paste0(oDir,"UCLex_",release,"_all_samples_variant_summary"),col.names=F,row.names=F,quote=F,sep="\t",append=T) 

		write.table(map, oMap, col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t", append=TRUE)  
		write.table(data.frame(cbind(map, annotations.snpStats$Obs,annotations.snpStats$Ref ) ) , oBim, col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t", append=TRUE) 
	}


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
