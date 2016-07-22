getArgs <- function() {
  myargs.list <- strsplit(grep("=",gsub("--","",commandArgs()),value=TRUE),"=")
  myargs <- lapply(myargs.list,function(x) x[2] )
  names(myargs) <- lapply(myargs.list,function(x) x[1])
  return (myargs)
}



myArgs <- getArgs()
print(myArgs)
if ('rootODir' %in% names(myArgs))  rootODir <- myArgs[[ "rootODir" ]]
if ('release' %in% names(myArgs))  release <- myArgs[[ "release" ]]

#######################################

oDir <- paste0(rootODir, "/UCLex_", release, "/")
print(oDir)
extCtrl.var <- read.table( paste0(oDir, "/Ext_ctrl_variant_summary") , header=T) 

hwe <- read.table(paste0(oDir, "gstats.hwe" ) , header=T ) 
hwe.pval <- .0001
hwe2 <- subset(hwe$SNP, hwe$TEST == "ALL" & hwe$P >= hwe.pval ) 

## some parameters
missingness.threshold <- 0 # snps higher than this are kept. 
min.maf <- 0.01
max.maf <- 0.5

clean.variants <- subset(extCtrl.var, extCtrl.var$Call.rate >= missingness.threshold) 
clean.variants$SNP <- clean.variants[,1]
percent.removed <- paste0("(", round(( nrow(extCtrl.var) - nrow(clean.variants)) / nrow(extCtrl.var)*100), "%)") 
message(paste(nrow(extCtrl.var) - nrow(clean.variants) , percent.removed, "variants are removed because of call rate") ) 
clean.variants <- clean.variants[clean.variants$SNP %in% hwe2, ]

annotations <- read.csv(file = paste0(oDir, "annotations.snpStat"), header=T, sep="\t") 
pass.snps <- annotations$clean.signature[which(annotations$FILTER == "PASS") ] 
qual<-data.frame(snp=rownames(annotations),qual=annotations$FILTER)
clean.out <- clean.variants$SNP[clean.variants$SNP %in% pass.snps]
write.table(qual, file = paste0(oDir, "snp_quality_score"), col.names=F, row.names=F, quote=F, sep="\t") 
write.table(clean.out, file = paste0(oDir, "Clean_variants"), col.names=F, row.names=F, quote=F, sep="\t") 


func <-  c("nonsynonymous SNV", "stopgain SNV", "nonframeshift insertion", "nonframeshift deletion", "frameshift deletion", "frameshift substitution", "frameshift insertion",  "nonframeshift substitution", "stoploss SNV")
lof <-  c("frameshift deletion", "frameshift substitution", "frameshift insertion",  "stoploss SNV")

funky <- annotations$clean.signature[annotations$ExonicFunc %in% unlist(func) ] 
rare <- subset(annotations$clean.signature, annotations$ESP6500si_ALL >= min.maf & annotations$ESP6500si_ALL <= max.maf & annotations$X1000g2012apr_ALL >= min.maf & annotations$X1000g2012apr_ALL <= max.maf ) 

#save(extCtrl.var, funky, rare, clean.variants, file = "tmp.RData")
funky.rare <- funky[funky %in% rare]
clean.variants.rare <- subset(extCtrl.var[,1] , extCtrl.var$Call.rate >= missingness.threshold & extCtrl.var$MAF >= min.maf & extCtrl.var$MAF <= max.maf) 

clean.funky <- clean.variants.rare[clean.variants.rare %in% funky.rare]
message(nrow(clean.funky))
write.table(clean.funky, file = paste0(oDir, "Clean_variants_func_rare" ), col.names=F, row.names=F, quote=F, sep="\t") 