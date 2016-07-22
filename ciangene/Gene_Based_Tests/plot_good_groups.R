getArgs <- function() {
  myargs.list <- strsplit(grep("=",gsub("--","",commandArgs()),value=TRUE),"=")
  myargs <- lapply(myargs.list,function(x) x[2] )
  names(myargs) <- lapply(myargs.list,function(x) x[1])
  return (myargs)
}

release <- 'February2015'
rootODir<-'/scratch2/vyp-scratch2/cian'
myArgs <- getArgs()

if ('rootODir' %in% names(myArgs))  rootODir <- myArgs[[ "rootODir" ]]
if ('release' %in% names(myArgs))  release <- myArgs[[ "release" ]]

######################
oDir <- paste0("/scratch2/vyp-scratch2/cian/UCLex_", release, "/LDAK_gene_tests_all_phenos/") 
group <- list.files(oDir, pattern="Lambiase_.0000001_.00000001", full.names=T) 
source("annotate_qqplot.R")


tech <- read.table(paste0(group[grep("with_kin_", group)], "/regressALL") , header=T)
genes<-tech$Gene_name[order(tech$LRT_P_Perm)][1:4]

pdf("tst.pdf")
par(mfrow=c(2,2))
for(i in 1:length(group))
{
	file <- read.table(paste0(group[i],"/regressALL"), header=T)
	tst <- qq.chisq(-2*log(file$LRT_P_Perm), df=2, x.max=30, main=basename(group)[i],cex.main=0.7)
	labelQQplot(file$Gene_name, tst, file$LRT_P_Perm) 
}
for(i in 1:length(group))
{
	file <- read.table(paste0(group[i],"/regressALL"), header=T)
	tst <- qq.chisq(-2*log(file$LRT_P_Perm), df=2, x.max=30, main=basename(group)[i],cex.main=0.7)
	labelQQplot(file$Gene_name, tst, file$LRT_P_Perm,GenesToLabel=genes) 
}
dev.off()
