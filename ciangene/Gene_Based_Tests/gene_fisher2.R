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
library(data.table)
oDir <- paste0(rootODir, "/UCLex_", release, "/")

library(SKAT)
library(MultiPhen)
library(parallel)

oFile <- "SNP.data"
if(!file.exists(oFile))
{
	dat <- read.table(paste0("/scratch2/vyp-scratch2/cian/UCLex_",release,"/annotations.snpStat"), header=T, sep="\t")
	maf.threshold <- 0.01
	dat$ESP6500si_ALL[is.na(dat$ESP6500si_ALL)] <- 0 
	dat$X1000g2012apr_ALL[is.na(dat$X1000g2012apr_ALL)] <- 0 
	dat <- subset(dat, dat$ESP6500si_ALL <= maf.threshold | dat$X1000g2012apr_ALL <= maf.threshold) 
	dat <- subset(dat, dat$FILTER == "PASS") 

	ext.ctrls <- read.table(paste0("/cluster/project8/vyp/cian/data/UCLex/UCLex_October2014/Single_variant/Clean_up/ExtCtrlMAFs", header=T, sep="\t") 
	ext.ctrls.filt <- subset(ext.ctrls, ext.ctrls$ExtCtrl_MAF <= maf.threshold) 
	
	dat <- dat[rownames(dat) %in% unlist(ext.ctrls.filt$SNP) , ] 

	snps <- data.frame(rownames(dat), dat$Gene, dat$ExonicFunc) 
	colnames(snps) <- c("SNPs", "Gene", "ExonicFunc") 
	write.table(snps, oFile, col.names=T, row.names=F, quote=F, sep="\t") 
}	else snps <- read.table(oFile, header=T, sep="\t") 



classes <- c("nonsynonymous SNV", "synonymous SNV", "stopgain SNV", "nonframeshift insertion", "nonframeshift deletion", "frameshift deletion", "frameshift substitution", "frameshift insertion", "unknown", "nonframeshift substitution stoploss SNV") 

func <-  c("nonsynonymous SNV", "stopgain SNV", "nonframeshift insertion", "nonframeshift deletion", "frameshift deletion", "frameshift substitution", "frameshift insertion",  "nonframeshift substitution", "stoploss SNV")
lof <-  c("frameshift deletion", "frameshift substitution", "frameshift insertion",  "stoploss SNV")  



spl <- t(data.frame(strsplit(as.character(snps$SNPs), "_"))) 

snps$Gene <- substr(snps$Gene, 1, 15 )
snps$Chr <- gsub(snps$SNPs, pattern = "_.*", replacement = "" )
snps$Chr <- gsub(snps$Chr, pattern = "X", replacement = 23 )
snps$BP <- as.numeric(spl[,2]) 

for (chr in 1:23)
{
	snps.small <- subset(snps, snps$Chr == chr)
	# tra <- data.frame(aggregate(snps$BP, list(snps$Chr) , summary) )
	bed <- data.frame(snps.small$Chr[1], min(snps.small$BP), max(snps.small$BP), snps.small$Chr[1] ) 
	write.table(bed, file = paste0(chr, "_ranges"), col.names=F, row.names=F, quote=F, sep="\t")
}


snps.func <- snps[ snps$ExonicFunc %in% unlist(func) , ]
write.table(snps.func, "SNPs.func", col.names=T, row.names=F, quote=F, sep="\t")
snps.lof <- snps[ snps$ExonicFunc %in% unlist(lof) , ]
write.table(snps.lof, "SNPs.lof", col.names=T, row.names=F, quote=F, sep="\t")



doFISHER <- function(data, t,  pheno, oBase, type)
{
	iFile <- paste0(data, t)
	if(file.exists(paste0(iFile, ".bed") ))
	{
		plink <- read.plink(iFile)

		majs <- apply(plink, 2, function(x) which(x == 0 ))
		mins <- apply(plink, 2, function(x) which(x == 2 ))

		for(i in 1:length(majs))  ## major alleles sometimes flipped, so this will flip em back. 
		{
			if( length(majs[[i]]) < length(mins[[i]]) ) 
			{
				plink[majs[[i]] ,i] <- 2
				plink[mins[[i]] ,i] <- 0
			}
		}

		bim <- read.table(paste0(iFile, ".bim"), header=F, sep="\t") 
		colnames(plink) <- bim[,2] 

		pheno <- read.table(pheno, header=F ,sep="\t")
		#pheno[grep("-9", pheno) ] <- NA
		cases <- pheno[pheno[,3] == 1 , 1]
		ctrls <- pheno[pheno[,3] == 0 , 1]

		type <- subset(type, type$Chr == bim[1,1])
		genes <- unique(type$Gene) 

		colNamesDat <- c("Gene", "FisherPvalue", "nb.mutations.cases", "nb.mutations.ctrls", "nb.cases", "nb.ctrls" , "Nb.variants")
		dat <- data.frame(matrix(nrow = length(genes), ncol = length(colNamesDat )  ) ) 
		colnames(dat)  <- colNamesDat
		dat[,1] <- genes

		for(i in 1:length(genes))
		{
			hit <- which(type$Gene == genes[i] )
			# hit <- unique(hit, grep( genes[i] , func.genes$Gene))
			hit.snps <- as.character(unlist(type$SNPs[hit] )) 
			snps <- as.matrix(plink[ ,  colnames(plink)  %in% hit.snps  ] ) 

			case.calls <- snps[ rownames(snps) %in% unlist(cases), ]
			ctrl.calls <- snps[ rownames(snps) %in% unlist(ctrls), ]

			number_mutations_cases <- sum( case.calls , na.rm=T )
			number_mutations_ctrls <- sum( ctrl.calls , na.rm=T ) 

			nb.cases <-  length(which(!is.na( unlist(case.calls)  )) ) 
			nb.ctrls <-  length(which(!is.na( unlist(ctrl.calls) )) ) 

			mean_number_case_chromosomes <- nb.cases * 2
			mean_number_ctrl_chromosomes <- nb.ctrls * 2

			nb.variants <- ncol(case.calls); message(nb.variants) 	; if(!is.integer(nb.variants)) nb.variants <- NA

			if (!is.na(number_mutations_cases)  & !is.na(number_mutations_ctrls)  )
			{
				if (nb.cases > 0 & nb.ctrls > 0)
				{
					fishertest <-  fisher.test((matrix(c(number_mutations_cases, 
									mean_number_case_chromosomes - number_mutations_cases, number_mutations_ctrls, 
									mean_number_ctrl_chromosomes - number_mutations_ctrls),
					                nrow = 2, ncol = 2)))
					dat[i,2] <- fishertest$p.value 
				}
			}
			dat$nb.mutations.cases[i] <- number_mutations_cases
			dat$nb.mutations.ctrls[i] <- number_mutations_ctrls
			dat$nb.cases[i] <- nb.cases 
			dat$nb.ctrls[i] <- nb.ctrls 
			dat$Nb.variants[i] <- nb.variants
		}

		oFile <- paste0(oBase, t)
		write.table(dat , oFile,  col.names=T, row.names=F, quote=F, sep="\t")
	}
}


# doFISHER <- function(data, t,  pheno, oBase, type)

ivfPheno <- '/cluster/project8/vyp/cian/data/UCLex/UCLex_October2014/Lambiase_case_control/support/IVF.pheno'
Lambiase_vs_UCLecx_Pheno <- "/cluster/project8/vyp/cian/data/UCLex/UCLex_October2014/Lambiase_case_control/support/Lambiase_vs_UCLex_phenotype_file"

data <- "/scratch2/vyp-scratch2/cian/UCLex_October2014/Lambiase_case_control/UCLex_"
#mclapply(1:23, mc.cores = 5, function(t)
for(t in 1:23)
{
	doFISHER (data , t, ivfPheno , "IVFvsUCLgene_func_FISHER_"  , snps.func ) 
	doFISHER( data , t, ivfPheno , "IVFvsUCLgene_lof_FISHER_", snps.lof) 

	doFISHER (data , t, Lambiase_vs_UCLecx_Pheno , "LambiaseVsUCLgene_func_FISHER_"  , snps.func ) 
	doFISHER( data , t, Lambiase_vs_UCLecx_Pheno , "LambiaseVsUCLgene_lof_FISHER_", snps.lof) 

}

#)




