getArgs <- function() {
  myargs.list <- strsplit(grep("=",gsub("--","",commandArgs()),value=TRUE),"=")
  myargs <- lapply(myargs.list,function(x) x[2] )
  names(myargs) <- lapply(myargs.list,function(x) x[1])
  return (myargs)
}

release <- 'July2015'
rootODir<-'/scratch2/vyp-scratch2/cian'

myArgs <- getArgs()

if ('rootODir' %in% names(myArgs))  rootODir <- myArgs[[ "rootODir" ]]
if ('release' %in% names(myArgs))  release <- myArgs[[ "release" ]]

######################
library(snpStats)
library(biomaRt) 

bDir <- paste0("/scratch2/vyp-scratch2/cian/UCLex_", release, "/")
minMaf <- c(0, 0.00001,0.01,0.1,0.2)
count.thresholds <- c(0,100,500,1000)
#count.thresholds<-c(0,0.00001,0.001,0.1)

iDir <- paste0(bDir, "FastLMM_Single_Variant_all_phenos/")

models<-c("base","tk_depth","perm")
AllFiles<-list.files(iDir,full.names=TRUE)
if(length(grep('sh',AllFiles))>0)AllFiles<- AllFiles[-grep('sh',AllFiles)]
if(length(grep('merged',AllFiles))>0)AllFiles<- AllFiles[-grep('merged',AllFiles)]
if(length(grep('pdf',AllFiles))>0)AllFiles<- AllFiles[-grep('pdf',AllFiles)]

Collate<-function(group,model)
{
	model.files<-AllFiles[grep(model,AllFiles)]
	files<-model.files[grep(paste0(group,"_"),model.files)]
	if(length(files)>0)
	{
		for(i in 1:length(files))
		{
			oFile<-paste0(iDir,model,"_",group,"_merged")
			#if(!file.exists(oFile))
			#{
				message(paste('Merging into',oFile))
				file<-read.table(files[i],header=T,sep="\t") 
				if(i==1)write.table(file,oFile,col.names=T, row.names=F, quote=F, sep="\t", append=F) 
				if(i>1)write.table(file,oFile,col.names=F, row.names=F, quote=F, sep="\t", append=T) 
			#} else message(paste(oFile,'already exists, so skipping'))
		}
	} else message("No files match input, wtf")
}

groups <- gsub(basename(AllFiles[grep("base",AllFiles)]), pattern = "_.*", replacement = "")
uniq.groups <- unique(groups)
start.group <- 1
end.group<-length(uniq.groups) 

prep<-TRUE
if(prep)
{
	for(i in start.group:end.group) 
	{
		for(mod in 1:length(models))
		{
			Collate(uniq.groups[i],models[mod])
		}
	}
}

files <- list.files(iDir, pattern = "merged", full.names=TRUE)

noKin <- files[grep('base',files)]
techKin <-files[grep('tk_depth',files)]
permy <- files[grep('perm',files)]

counts <- list.files(paste0(bDir, "Single_variant_tests/") , pattern = "assoc$", full.names=TRUE) 
noCovar <- list.files(paste0(bDir, "Single_variant_tests/") , pattern = ".*no.*adjusted", full.names=TRUE)
techCovar <- list.files(paste0(bDir, "Single_variant_tests/") , pattern = ".*tech.*adjusted", full.names=TRUE)

annotations <- read.csv(paste0(bDir, "annotations.snpStat"), header=TRUE, sep="\t")

ex.ctrl.dir<-paste0(bDir, "/External_Control_data/") 



oFile <- paste0(iDir, "SingleVariant_qqplots_counts.pdf")
pdf(oFile)
par(mfrow=c(2,2), cex.main=0.8)
for(i in 1:length(uniq.groups))
{
	iBase <- noKin[grep(paste0("_",uniq.groups[i],"_"), noKin) ]
	iTech <- techKin[grep(paste0("_",uniq.groups[i],"_"), techKin) ]
	#iCovar <- noCovar[grep(paste0(uniq.groups[i],"_"), noCovar) ]
	#iTk <- techCovar[grep(paste0(uniq.groups[i],"_"), techCovar) ]
	iPerm <- permy[grep(paste0(uniq.groups[i],"_"), permy) ]
#	iRes<-inRes[grep(paste0("_",uniq.groups[i],"_"),inRes) ]
	inputs <- c(iBase, iTech, iPerm)#,iRes)
	if(length(which(file.exists(inputs)))==3)
	{
		message("Now processing ", uniq.groups[i])
		base <- read.table(iBase, header=T, stringsAsFactors=F)
		tech <- read.table(iTech, header=T, stringsAsFactors=F)
		tech.small <- data.frame(SNP=tech$SNP, TechKinPvalue=tech$Pvalue)
		perm <- read.table(iPerm, header=T, stringsAsFactors=F)
		perm.small <- data.frame(SNP=perm$SNP, permPvalue=perm$Pvalue)
#		res <- read.table(iRes, header=T, stringsAsFactors=F)
#		res.small <- data.frame(SNP=res$SNP, resPvalue=res$Pvalue)

		results.merged2 <- merge(base, tech.small, by = "SNP",all=T)
#		results.merged2<-merge(results.merged, res.small, by = "SNP")
		results.merged3 <- merge(results.merged2,perm.small,by='SNP',all=T)
		results.merged.anno2 <- merge(results.merged3, annotations, by.x = "SNP", by.y ="clean.signature",all=T)

	#	noCov <- read.table(iCovar, header=T, stringsAsFactors=F)
	#	noCov.small<-data.frame(SNP=noCov$SNP, noCov$UNADJ)
	#	tk <- read.table(iTk, header=T, stringsAsFactors=F)
#		tk.small<-data.frame(SNP=tk$SNP,tk$UNADJ)
#		covariate.pvalues<-merge(noCov.small, tk.small,by='SNP',all=T)

	#	results.merged.anno2 <- merge(results.merged.anno, covariate.pvalues,by='SNP',all=T)

		ex.ctrl.file<-paste0(ex.ctrl.dir,uniq.groups[i],'_ex_ctrls.frq')
		system( paste('tr -s " " <',ex.ctrl.file , '>', paste0(ex.ctrl.file,'_clean') ) ) 
		extCtrl <- read.table(paste0(ex.ctrl.file,'_clean'),header=T,sep=" ") 
		extCtrl.small <- data.frame(extCtrl$SNP, extCtrl$MAF)
		colnames(extCtrl.small) <-  c("SNP", "ExtCtrl_MAF") 
		extCtrl.small$ExtCtrl_MAF[is.na(extCtrl.small$ExtCtrl_MAF)] <- 0 
		message('here')
		results.merged.anno.extCtrl <- merge(results.merged.anno2, extCtrl.small, by = "SNP",all=T)
		#results.merged.anno.extCtrl$tk.UNADJ[which(is.infinite(results.merged.anno.extCtrl$tk.UNADJ) )] <- 1
		hit <- counts[grep(uniq.groups[i], counts)][1]
		if(file.exists(hit))
		{
		current.counts <- read.table(hit, header=T) 
		final <- merge(results.merged.anno.extCtrl, current.counts, by = "SNP",all=T)
		final<-final[order(as.numeric(as.character(final$TechKinPvalue))),]
		oFile<-paste0(iDir, uniq.groups[i], "_final");message(paste("Writing to", oFile)) 

		ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
		filter="ensembl_gene_id"
		attributes =  c("ensembl_gene_id", "external_gene_name",  "phenotype_description")
		print("Done tests, now getting gene names and phenotype")
		dat.anno <- getBM(attributes= attributes , filters = filter , values = final$Gene , mart = ensembl)
		skat.annotated <- dat.anno[!duplicated(dat.anno$external_gene_name),]
		merged<-merge(final,skat.annotated,by.x="Gene",by.y="ensembl_gene_id",all.x=T)
		print(head(merged))
		merged<-merged[order(merged$TechKinPvalue),]
		write.table(merged, oFile, col.names=T, row.names=F, quote=F, sep="\t")
	#	lapply(minMaf, function(x)
		message("Now plotting ", uniq.groups[i])
		lapply(count.thresholds,function(x)
		{
	#		dat <- data.frame(subset(final, final$ExtCtrl_MAF >= x )) 
			dat <- subset(final, final$C_U >= x )
			message(nrow(dat),' rows with counts greater than ',x)
			if(nrow(dat) > 0 )
			{
				qq.chisq(-2*log(as.numeric(dat$Pvalue)),df=2,x.max=30,pvals=T, main = paste(uniq.groups[i], x, "noKin",nrow(dat)))
				qq.chisq(-2*log(as.numeric(dat$TechKinPvalue)),df=2,x.max=30, pvals=T, main = paste(uniq.groups[i], x, "TechKin", nrow(dat)))
			}   
		}
		)
	} ## file.xists(hit)
} else message("Skipping " , uniq.groups[i]) # file.exists(inputs)
}
dev.off()
source("summarise_single_variant_tests.R") 
