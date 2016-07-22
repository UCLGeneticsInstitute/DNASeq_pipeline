gene.list <- "/SAN/biomed/biomed14/vyp-scratch/cian/LDAK/genesldak_ref.txt" 
gene.list <- read.table(gene.list, header=F, sep="\t" ) 

release <- "August"

dir <- paste0("/cluster/project8/vyp/cian/data/UCLex/UCLex_", release, "/Scripts/Gene_based_tests/") 
dirs <- list.files(dir , pattern ="_kin_m") 

prep <- FALSE	
if(prep){
	for(i in 1:length(dirs)) { 
		system( paste("rm",  paste0(dirs[i], "/effects_ALL") ) )
		system( paste("cat", paste0(dirs[i], "/effects*") , ">> ", paste0(dirs[i], "/effects_ALL") ) )
	} 
} 


getWeights <- function(gene, cohort) { 
	dir <- paste0("/cluster/project8/vyp/cian/data/UCLex/UCLex_", release, "/Scripts/Gene_based_tests/") 
	dirs <- list.files(dir , pattern = paste0("_", cohort, "_")) 

	files <- list.files(dirs, pattern = "effects_ALL", full.names=T) 
	names <- gsub(files, pattern = "_missingness_0.9", replacement ="")

	for(i in 1:length(files)) { 
	file <- read.table(files[i], header=T, sep=" ", stringsAsFactors=F) 
	genes <- read.table(paste0(dirs[i], "/gene_details.txt") , header=T, sep=" ", skip = 2) 
	gene_nb <- genes$Gene_Number[ grepl(gene, genes$Gene_Name) ] 
	file.small <- subset(file, file$Gene_Number == gene_nb) 

	if(i ==1) { 
		dat <- data.frame(matrix(nrow=nrow(file.small),ncol=length(files)+2) ) 
		dat[,1] <- gene
		dat[,2] <- file.small$Predictor
	} 
	dat[,i+2] <- file.small$Effect	
	
	} 
names(dat) <- c("Gene", "Predictor", unlist( gsub(names, pattern = "/.*", replacement = "")  ) )
return(dat)
} 


dat <- getWeights("SAMD11", "Hardcastle")
