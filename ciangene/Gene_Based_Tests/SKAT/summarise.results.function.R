suppressPackageStartupMessages(library(snpStats))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(xtable))
suppressPackageStartupMessages(library(data.table))

ldak='/cluster/project8/vyp/cian/support/ldak/ldak'
source('/cluster/project8/vyp/cian/scripts/r/annotate_qqplot.R') ## Modified qqplot so I can label specific genes. 

lof <-  c("frameshift deletion", "frameshift substitution", "frameshift insertion",  "stoploss SNV", "splicing"
		,"stopgain SNV","exonic;splicing"
		)
snps<-read.table("/SAN/vyplab/UCLex/mainset_July2016/cian/Annotations/func.tab",sep='\t',header=T) 
snps.filt<-snps[snps$ExonicFunc %in% lof & snps$ExAC_MAF < 0.0001,]
snps.filt<-snps.filt[!duplicated(snps.filt$SNP),]

data<-'/SAN/vyplab/UCLex/mainset_July2016/cian/allChr_snpStats_out'


summarise<-function(dir,genes=NULL,outputDirectory='Results',plot=TRUE,Title=basename(outputDirectory))
{
	outputDirectory<-paste0(outputDirectory,'/')
	if(!file.exists(outputDirectory)) dir.create(outputDirectory)

	file.copy(paste0(dir,'qc/case_pca.plot.pdf'),paste0(outputDirectory,'Case_PCA_plots.pdf')) 
	message("Added PCA plot of cases to outputDirectory")
	file.copy('/SAN/vyplab/UCLex/support/SKAT/Explanation_SKAT_CSV_Columns.docx',paste0(outputDirectory,'Explanation_SKAT_CSV_Columns.docx')) 

	files<-list.files(dir,pattern='csv',full.names=T,recursive=T)
	files<-files[grep("chr",files)]

	oFile<-paste0(outputDirectory,Title,'_SKAT_results.csv') 
	if(file.exists(oFile))file.remove(oFile)
	message("Merging results found in:")
	for(i in 1:length(files))
	{
		file<-read.csv(files[i],header=T)
		names<-colnames(file)
		#oFile<-paste0(dataDir,files[i,2],'_macular_dystrophy_results.csv') 
		print(oFile)
		write.table(file,oFile,col.names=!file.exists(oFile),row.names=F,quote=T,sep=',',append=T)
	}

	file<-read.csv(oFile,header=F)
	file.remove(oFile)
	colnames(file)<-names
	file$SKATO<-as.numeric(as.character(file$SKATO))
	file<-file[order(file$SKATO),]
	file<-file[-grep("nb.cases",file$nb.cases),]
	write.table(file,oFile,col.names=T,row.names=F,quote=T,sep=',')

	MergeFiles<-function(files,out)
	{
		if(file.exists(out))file.remove(out)
		for(carry in 1:length(files))
		{
			c.file<-unique(read.csv(files[carry],header=F,sep='\t'))
			write.table(c.file, out,col.names=F,row.names=F,quote=T,sep=',',append=T)
		}
		outFile<-read.csv(out,header=FALSE)
		return(outFile)
	}
	
	carriers<-list.files(dir,pattern='case_carriers',full.names=T,recursive=T)
	carry.oFile<-paste0(outputDirectory,Title,'_case_carriers.csv') 
	carriers<-MergeFiles(carriers,carry.oFile)
	
	ctrl.carriers<-list.files(dir,pattern='ctrl_carriers',full.names=T,recursive=T)
	ctrl.carry.oFile<-paste0(outputDirectory,Title,'_control_carriers.csv') 
	ctrl.carriers<-MergeFiles(ctrl.carriers,ctrl.carry.oFile)

	by.snp<-list.files(dir,pattern='SKAT_results_by_SNP',full.names=T,recursive=T)
	by.snp.out<-paste0(outputDirectory,Title,'_results_by_SNP_cases.csv') 
	by.snp<-MergeFiles(by.snp,by.snp.out)
	colnames(by.snp)[1]<-'SNP'

	CompoundHets_cases<-list.files(dir,pattern='CompoundHets_cases',full.names=T,recursive=T)
	CompoundHets_cases.oFile<-paste0(outputDirectory,Title,'_CompoundHets_cases.csv') 
	CompoundHets_cases<-MergeFiles(CompoundHets_cases,CompoundHets_cases.oFile)
	
	CompoundHets_ctrls<-list.files(dir,pattern='CompoundHets_ctrls',full.names=T,recursive=T)
	CompoundHets_ctrls.oFile<-paste0(outputDirectory,Title,'_CompoundHets_ctrls.csv') 
	CompoundHets_ctrls<-MergeFiles(CompoundHets_ctrls,CompoundHets_ctrls.oFile)

	SummariseCpdHets<-function(file)
	{
		nb.snp.rows<-ncol(file)-5
		file$LOF.rare.snps<-0 

		for(row in 1:nrow(file))
		{	
			snps<-unlist(file[row,3:(3+nb.snp.rows)]) 
			nb.snps<-length(which(snps %in% snps.filt$SNP))
			file$LOF.rare.snps[row]<-nb.snps
		}
		file<-subset(file,file$LOF.rare.snps>1)
		return(file) 
	}

	cpd.hets.cases<-SummariseCpdHets(CompoundHets_cases)
	cpd.hets.ctrls<-SummariseCpdHets(CompoundHets_ctrls)
	if(nrow(cpd.hets.cases)>0)write.table(cpd.hets.cases,CompoundHets_cases.oFile,col.names=T,row.names=F,quote=F,sep='\t') else file.remove (CompoundHets_cases.oFile)
	if(nrow(cpd.hets.ctrls)>0)write.table(cpd.hets.ctrls,CompoundHets_ctrls.oFile,col.names=T,row.names=F,quote=F,sep='\t') else file.remove (CompoundHets_ctrls.oFile)

	if(plot)
	{

		qqplot.name<-paste0(paste0(outputDirectory,'qqplots.pdf')) 
		message(paste('QQplot stored in:',qqplot.name))
		pdf(qqplot.name)


		file<-file[which(file$SKATO>0),]

		if(!is.null(genes))
		{
			if( file.exists(file.path(genes))) genes<-read.table(genes,header=FALSE)[,1]
			file$Candidate<-FALSE
			file$Candidate[file$Symbol %in% genes]<-TRUE
			file<-data.frame(cbind(file[,ncol(file)],file[,1:(ncol(file)-1)]))
			colnames(file)[1]<-'Candidate'
			tst <- qq.chisq(-2*log(file$SKATO), df=2, x.max=50,pvals=TRUE,main=Title)
			labelQQplot(pval.labels=file$Symbol,qqplotOut=tst,pvals.raw=file$SKATO, genelist= genes) 
		} else{
			sig.genes<-subset(gene.pvals,gene.pvals$SKATO <0.0001)$Symbol
			tst <- qq.chisq(-2*log(file$SKATO), df=2, x.max=50,pvals=TRUE,main=Title)
			labelQQplot(pval.labels=file$Symbol,qqplotOut=tst,pvals.raw=file$SKATO, genelist= sig.genes )
		}
		dev.off()
	}

	file$Carriers<-0
	for(snp in 1:nrow(carriers))
	{
		hit<-grep(carriers$V1[snp],file$SNPs)
		if(length(hit)>0)
		{
			for(snp.hit in 1:length(hit))
			{
				if(file$Carriers[hit[snp.hit]]!=0)file$Carriers[hit[snp.hit]]<-paste0(file$Carriers[hit[snp.hit]],';',carriers$V2[snp])
				if(file$Carriers[hit[snp.hit]]==0)file$Carriers[hit[snp.hit]]<-carriers$V2[snp]
			}
		}
	}
	write.table(file,paste0(outputDirectory,Title,'_SKAT_processed.csv'),col.names=T,row.names=F,quote=T,sep=',',append=F)
	
	filt<-file[1:10,]
	rownames(filt)<-1:nrow(filt)

	filt.xtable<-xtable(filt,caption=paste(Title,"SKAT top genes") ,digits=2, display = c(rep("s",3),rep("E",26),rep("s",3)))
	htmlOut<-paste0(outputDirectory,Title,"_SKAT.html")
	print(htmlOut)
	print.xtable(filt.xtable, type="html",file=htmlOut,scalebox=.7)

	tt<-merge(by.snp[,1:10],snps.filt,by='SNP')
	if(nrow(tt)>0)
	{
		tt$Carriers<-NA
		for(i in 1:nrow(tt))
		{
			hit<-grep(tt$SNP[i],carriers[,1])
			tt$Carriers[i]<-paste0(carriers[hit,2],collapse=';')
		}
		case.ctrl.carriers<-subset(tt,tt$case.snp.hets>1 | tt$case.snp.homs > 1 )
		#if(nrow(case.ctrl.carriers)>0)write.table(case.ctrl.carriers,paste0(outputDirectory,Title,'_BiallelicSNPs.csv'),col.names=T,row.names=F,quote=T,sep=',',append=F)
		#if(nrow(case.ctrl.carriers)==0)write.table('None Found',paste0(outputDirectory,Title,'_BiallelicSNPs.csv'),col.names=T,row.names=F,quote=T,sep=',',append=F)
	}
} # summarise