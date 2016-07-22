CustomPheno<-function(basePheno,remove)
{
	# basePheno is initial pheno file, eg the one from make_phenotype_file.R
	# remove is a string
	for(i in 1:length(remove))
	{
		if(i==1)filtPheno<-basePheno
		hit<-grep(remove[i],basePheno[,1])
		if(length(hit)==0){ warning(paste("No samples called" ,remove[i],"found"))}else{
			filtPheno[hit,3:ncol(filtPheno)]<-NA
		}
	}
return(filtPheno)
}

syrris<-CustomPheno(pheno,remove=c("Lambiase","test")  ) 
lambiase<-CustomPheno(pheno,remove=c("Syrris","test")  ) 