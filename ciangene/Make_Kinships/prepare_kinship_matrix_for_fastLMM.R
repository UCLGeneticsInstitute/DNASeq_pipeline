getArgs <- function() {
  myargs.list <- strsplit(grep("=",gsub("--","",commandArgs()),value=TRUE),"=")
  myargs <- lapply(myargs.list,function(x) x[2] )
  names(myargs) <- lapply(myargs.list,function(x) x[1])
  return (myargs)
}

release <- 'February2015'

myArgs <- getArgs()

if ('rootODir' %in% names(myArgs))  rootODir <- myArgs[[ "rootODir" ]] else rootODir <- paste0("/scratch2/vyp-scratch2/cian/UCLex_", release, "/") 
if ('release' %in% names(myArgs))  release <- myArgs[[ "release" ]]
#########################################

bDir <- paste0(rootODir, "/UCLex_", release, "/")

kin=as.matrix(read.table(paste0(bDir, "TechKin.grm.raw"))) 
fam=as.matrix(read.table(paste0(bDir, "allChr_snpStats.fam")))
fam=as.matrix(read.table(paste0(bDir, "TechKin.grm.id")))

names=paste(fam[,1],fam[,2])
names2="var";for(i in
1:nrow(fam)){names2=c(names2,paste(fam[i,1],fam[i,2]))}

write.table(rbind(names2), file = paste0(bDir, "TechnicalKinship_Fastlmm") ,row=F,col=F,quote=F,sep="\t")
write.table(cbind(names,kin),file = paste0(bDir, "TechnicalKinship_Fastlmm"),row=F,col=F,quote=F,sep="\t",append=T)

