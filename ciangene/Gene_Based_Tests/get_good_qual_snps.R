oFile <- "SNP.data"
if(!file.exists(oFile))
{
	dat <- read.table("/scratch2/vyp-scratch2/cian/UCLex_October2014/annotations.snpStat", header=T, sep="\t")
	maf.threshold <- 0.01
	dat$ESP6500si_ALL[is.na(dat$ESP6500si_ALL)] <- 0 
	dat$X1000g2012apr_ALL[is.na(dat$X1000g2012apr_ALL)] <- 0 
	dat <- subset(dat, dat$ESP6500si_ALL >= maf.threshold | dat$X1000g2012apr_ALL >= maf.threshold) 
	dat <- subset(dat, dat$FILTER == "PASS") 

	ext.ctrls <- read.table("/cluster/project8/vyp/cian/data/UCLex/UCLex_October2014/Single_variant/Clean_up/ExtCtrlMAFs", header=T, sep="\t") 
	ext.ctrls.filt <- subset(ext.ctrls, ext.ctrls$ExtCtrl_MAF >= maf.threshold) 
	
	dat <- dat[rownames(dat) %in% unlist(ext.ctrls.filt$SNP) , ] 

	snps <- data.frame(rownames(dat), dat$Gene, dat$ExonicFunc) 
	colnames(snps) <- c("SNPs", "Gene", "ExonicFunc") 
	write.table(snps, oFile, col.names=T, row.names=F, quote=F, sep="\t") 
}	else snps <- read.table(oFile, header=T, sep="\t") 

classes <- c("nonsynonymous SNV", "synonymous SNV", "stopgain SNV", "nonframeshift insertion", "nonframeshift deletion", "frameshift deletion", "frameshift substitution", "frameshift insertion", "unknown", "nonframeshift substitution stoploss SNV") 

func <-  c("nonsynonymous SNV", "stopgain SNV", "nonframeshift insertion", "nonframeshift deletion", "frameshift deletion", "frameshift substitution", "frameshift insertion",  "nonframeshift substitution", "stoploss SNV")
lof <-  c("frameshift deletion", "frameshift substitution", "frameshift insertion",  "stoploss SNV")  


snps.func <- snps[ snps$ExonicFunc %in% unlist(func) , ]
write.table(snps.func, "SNPs.func", col.names=T, row.names=F, quote=F, sep="\t")
snps.lof <- snps[ snps$ExonicFunc %in% unlist(lof) , ]
write.table(snps.lof, "SNPs.lof", col.names=T, row.names=F, quote=F, sep="\t")

