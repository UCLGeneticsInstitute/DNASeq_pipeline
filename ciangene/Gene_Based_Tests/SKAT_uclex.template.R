library(SKAT)
library(snpStats)
library(HardyWeinberg)

oFile<-paste0(oDir,pheno.matching[phen,1],'_snps')
if(file.exists(oFile))file.remove(oFile)

carrierFile<-paste0(oDir,pheno.matching[phen,1],'_snps_carriers')
if(file.exists(carrierFile))file.remove(carrierFile)

results$case.count<-0
results$ctrl.count<-0
results$case.maf<-0
results$ctrl.maf<-0
results$total.maf<-0
results$nb.case.homs<-0
results$nb.case.hets<-0
results$nb.ctrl.homs<-0
results$nb.ctrl.hets<-0 
results$Chr<-0
results$Start<-0
results$End<-0 


write.table( data.frame( colnames(clean.pheno.snps), current.pheno),paste0(oDir,pheno.matching[phen,1],'_case_control_list'), col.names=F,row.names=F,quote=F,sep='\t') 


for(gene in 1:nb.genes)
{
	gene.snps<-good.genes.data$SNP[ grep(uniq.genes[gene],good.genes.data$ENSEMBL) ]

	gene.data<- data.frame(t(data.frame(strsplit(gene.snps,'_')))) 
	gene.chr<-as.numeric(unique(gene.data[,1]))
	gene.start<-min(as.numeric(unique(gene.data[,2])))
	gene.end<-max(as.numeric(unique(gene.data[,2])))

	gene.snp.data<-clean.pheno.snps[ rownames(clean.pheno.snps) %in% gene.snps ,]
	nb.snps.in.gene<-nrow(gene.snp.data)
	print(paste(nb.snps.in.gene,'snps in', uniq.genes[gene]))

	test.gene<-unique(snp.gene$ENSEMBL[snp.gene$SNP %in% rownames(gene.snp.data)])  ## match to uniq. genes as a check, 

	if(length(unique(test.gene))>1) stop ("SNPs span multiple genes")
	if(test.gene!=uniq.genes[gene]) stop ("Genes not sorted correctly")
	if(results$ENSEMBL[gene]!=uniq.genes[gene]) stop ("Genes not sorted correctly")

	results$Chr[gene]<-gene.chr
	results$Start[gene]<-gene.start
	results$End[gene]<-gene.end

	if(nb.snps.in.gene>0)
	{
		case.snps<-gene.snp.data[,current.pheno==1]
		ctrl.snps<-gene.snp.data[,current.pheno==0]
		results$case.count[gene]<-sum(case.snps,na.rm=T)
		results$ctrl.count[gene]<-sum(ctrl.snps,na.rm=T)

		##these counts are for gene total
		if(sum(case.snps,na.rm=T)>0)results$case.maf[gene]<-signif(maf(as.numeric(unlist(table(unlist(case.snps))))),2) 
		
		case.homs<-length(grep(2,case.snps))
		if(case.homs>0) results$nb.case.homs[gene]<-case.homs

		case.hets<-length(grep(1,case.snps))
		if(case.hets>0)results$nb.case.hets[gene]<-case.hets

		### now do ctrls	
		if(sum(ctrl.snps,na.rm=T)>0)results$ctrl.maf[gene]<-signif(maf(as.numeric(unlist(table(unlist(ctrl.snps))))),2) 

		ctrl.homs<-length(grep(2,ctrl.snps))
		if(ctrl.homs>0)results$nb.ctrl.homs[gene]<-ctrl.homs
	
		ctrl.hets<-length(grep(1,ctrl.snps)) 
		if(ctrl.hets>0) results$nb.ctrl.hets[gene]<-ctrl.hets

		if(sum(gene.snp.data,na.rm=T)>0)results$total.maf[gene]<-signif(maf(as.numeric(unlist(table(unlist(gene.snp.data))))),2) 

		obj<-SKAT_Null_Model(current.pheno ~ 1, out_type="D")
		#results$SKAT[gene] <- SKAT(t(as.matrix(gene.snp.data)) , obj, missing_cutoff=0.4, estimate_MAF=2)$p.value 
		results$SKATO[gene] <- SKAT(t(as.matrix(gene.snp.data)) , obj, missing_cutoff=0.8, estimate_MAF=2,method="optimal.adj")$p.value
		results$nb.snps[gene] <- nb.snps.in.gene


		## these counts are for each snp in gene separately
		case.snp.hets<-apply(case.snps,1,function(x) length(grep(1,x)))
		case.snp.homs<-apply(case.snps,1,function(x) length(grep(2,x)))
		maf.snp.cases<-apply(case.snps,1, function(x) signif(maf(as.numeric(unlist(table(unlist(x))))),2)  )


		ctrl.snp.hets<-apply(ctrl.snps,1,function(x) length(grep(1,x)))
		ctrl.snp.homs<-apply(ctrl.snps,1,function(x) length(grep(2,x)))
		maf.snp.ctrls<-apply(ctrl.snps,1, function(x) signif(maf(as.numeric(unlist(table(unlist(x))))),2)  )

		snp.out<-data.frame( rownames(gene.snp.data),case.snp.hets,case.snp.homs, maf.snp.cases,ctrl.snp.hets,ctrl.snp.homs, maf.snp.ctrls,results[gene,1:8]) 
		write.table( snp.out, oFile, col.names=!file.exists(oFile),row.names=F,quote=F,sep='\t',append=T)
	
		case.carriers<-  rownames( data.frame(unlist(apply(case.snps,1,function(x) which(x>0 )))))
		ctrl.carriers<-  rownames( data.frame(unlist(apply(ctrl.snps,1,function(x) which(x>0 )))))
		carriers<- t(data.frame(c(case.carriers,ctrl.carriers) ) )
		write.table(carriers,carrierFile, col.names=FALSE,row.names=F,quote=F,sep='\t',append=T)

	}
}

results <- results[order(results[,2]), ]
results.out <- paste0(oDir,pheno.matching[phen,1],'_skat.csv')
qqplot.out <- paste0(rootODir,'/plots/',pheno.matching[phen,1],'_skat_QQplot.png')
write.table(results,results.out,col.names=T,row.names=F,quote=F,sep='\t')
png(qqplot.out)
qq.chisq(-2*log(results$SKATO), df=2, x.max=30, main=paste(pheno.matching[phen,1],'SKAT'),cex.main=0.7)	
dev.off() 