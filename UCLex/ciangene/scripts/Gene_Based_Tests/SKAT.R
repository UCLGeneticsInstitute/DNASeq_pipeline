library(SKAT)
library(MultiPhen)
library(parallel)

## Ran at bottom of script. give it a phenotype file, data and a list of snps to include and it'll SKAT them for you. 

oFile <- "SNP.data" ## This is where I filter the SNP list to start. I still need to include the case/ctrl call rate filter for gene based tests 
if(!file.exists(oFile))
{
	dat <- read.table("/scratch2/vyp-scratch2/cian/UCLex_October2014/annotations.snpStat", header=T, sep="\t")
	maf.threshold <- 0.01
	dat$ESP6500si_ALL[is.na(dat$ESP6500si_ALL)] <- 0 
	dat$X1000g2012apr_ALL[is.na(dat$X1000g2012apr_ALL)] <- 0 
	dat <- subset(dat, dat$ESP6500si_ALL <= maf.threshold | dat$X1000g2012apr_ALL <= maf.threshold) 
	dat <- subset(dat, dat$FILTER == "PASS") 

	ext.ctrls <- read.table("/cluster/project8/vyp/cian/data/UCLex/UCLex_October2014/Single_variant/Clean_up/ExtCtrlMAFs", header=T, sep="\t") 
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



doSKAT <- function(data, t,  pheno, oBase, type)
{
	iFile <- paste0(data, t)
	if(file.exists(paste0(iFile, ".bed") ))
	{
		plinky.func <- read.plink(iFile)
		
		bim <- read.table(paste0(iFile, ".bim"), header=F, sep="\t") 
		colnames(plinky.func) <- bim[,2] 

		pheno <- read.table(pheno, header=F ,sep="\t")#[,3]
	
		plinky.func <- plinky.func[rownames(plinky.func) %in% pheno[,1], ]

		pheno <- pheno[,3] 
		pheno <- pheno[pheno != '-9']

		type <- subset(type, type$Chr == bim[1,1])
		genes <- unique(type$Gene) 
		func <- data.frame(matrix(nrow= length(genes), ncol = 2))
		func[,1] <- as.character(genes) 


		for(i in 1:length(genes))
		{
			hit <- which(type$Gene == genes[i] )
			# hit <- unique(hit, grep( genes[i] , func.genes$Gene))
			hit.snps <- as.character(unlist(type$SNPs[hit] )) 
			snps <- as.matrix(plinky.func[ ,  colnames(plinky.func)  %in% hit.snps  ] ) 
		#	snps <- as.matrix(snps[!rownames(snps) %in% remove,])
			if(ncol(snps) > 0 )
			{
			obj<-SKAT_Null_Model(pheno[!is.na(pheno)] ~ 1, out_type="D")
			func[i,2] <- SKAT(snps, obj, missing_cutoff=0.4, estimate_MAF=2, kernel = "linear")$p.value
			message(func[i,])
			} 
		}

		func <- func[order(func[,2]), ]

		oFile <- paste0(oBase, t)
		write.table(func, oFile,  col.names=F, row.names=F, quote=F, sep="\t")
	}
}

# doFISHER <- function(data, t,  pheno, oBase, type)



ivfPheno <- '/cluster/project8/vyp/cian/data/UCLex/UCLex_October2014/Lambiase_case_control/support/IVF.pheno'
Lambiase_vs_UCLecx_Pheno <- "/cluster/project8/vyp/cian/data/UCLex/UCLex_October2014/Lambiase_case_control/support/Lambiase_vs_UCLex_phenotype_file"


data <- "/scratch2/vyp-scratch2/cian/UCLex_October2014/Lambiase_case_control/UCLex_"

for(t in 1:23)
{
	doSKAT (data , t, ivfPheno , "IVFvsUCLgene_func_SKAT_"  , snps.func ) 
	doSKAT( data , t, ivfPheno , "IVFvsUCLgene_lof_SKAT_", snps.lof) 

	doSKAT (data , t, Lambiase_vs_UCLecx_Pheno , "LambiaseVsUCLgene_func_SKAT_"  , snps.func ) 
	doSKAT( data , t, Lambiase_vs_UCLecx_Pheno , "LambiaseVsUCLgene_lof_SKAT_", snps.lof) 

}



