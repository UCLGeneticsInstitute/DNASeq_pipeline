oDir <- "Results"
if(!file.exists(oDir)) dir.create(oDir)

files <- list.files(, pattern = "_kin", full.names=T)
short <- gsub(files, pattern = ".*_kin", replacement="")
oFiles <- gsub(gsub(short, pattern = ".*pheno_", replacement = ""), pattern = "_.*", replacement = "") 


summariseModels <- function(Val, by.cohort=TRUE) { 
		for(i in 1:length(files)){
			target <- paste0(files[i], "/regressALL")
			if(file.exists(target)) { 
				file <- read.table(target, header=T, sep=" ")
				if(i==1) { 
				dat <- data.frame(matrix(nrow=nrow(file), ncol = length(files) + 1) ) 
				colnames(dat) <- c("Gene", unlist(basename(files))) 
				dat[,1] <- file$Gene_name
				}
				dat[,i+1] <- file[,which(colnames(file) == Val)]
			} 
		}
		
		if(by.cohort) { 
		cohort <- gsub(gsub(colnames(dat)[2:ncol(dat)], pattern = ".*pheno_", replacement = "") , pattern = "_.*", replacement ="") 
		cohort.uniq <- unique(cohort)
		listy <-  vector("list", length(cohort.uniq))
			for(t in 1:length(cohort.uniq)) listy[[t]] <-  data.frame(dat[,1], dat[ , which(cohort %in% cohort.uniq[t])+1 ] )
 		} 
	return(listy)	
} 


scoreP <- summariseModels("Score_P") ## feed this the name of the column you want. It will return a list, each item of the list a dataframe of this value by cohort per model.



plotPvalues <- function(data, cohort) { 
	for(i in 1:length(data)) { 
	hit <- grepl(cohort, colnames(data[[i]]) ) 
	if(length(which(hit)) > 0 ) break
	} 
	dat <- data[[i]]
	dat$base_res_diff <- dat[, paste0( "no_kin_maf_0.00001_pheno_", cohort, "_missingness_0.9") ] - dat[, paste0( "Tech_with_kin_maf_0.00001_pheno_", cohort, "_missingness_0.9") ]
	dat <- dat[order(dat$base_res_diff),]

	oFile <- paste0(cohort, ".pdf") 
	pdf(oFile)
	par(mfrow=c(2,2)) 
	hist(dat[, paste0( "no_kin_maf_0.00001_pheno_", cohort, "_missingness_0.9") ]   ,	
		 main = paste(cohort, "No kinships") , ylim = c(0,5000) , xlab = "Score test pvalues" )
	hist(dat[, paste0( "Depth_with_kin_maf_0.00001_pheno_", cohort, "_missingness_0.9") ]	,	
		 main = paste(cohort,"DepthKin") , ylim = c(0,5000) , xlab = "Score test pvalues" )
	hist(dat[, paste0( "Tech_with_kin_maf_0.00001_pheno_", cohort, "_missingness_0.9") ]	,	
		 main = paste(cohort,"TechKin") , ylim = c(0,5000)  , xlab = "Score test pvalues" )
	hist(dat[, paste0( "Tech_with_kin_maf_0.00001_pheno_" , cohort, "_missingness_0.9_res") ]	,
		 main = paste(cohort, "TechKin and DepthKin") , ylim = c(0,5000)  , xlab = "Score test pvalues" )
	dev.off()
} 

plotPvalues(scoreP, "Levine" ) ## simple histogram of some pvals
plotPvalues(scoreP, "Hardcastle" ) 











