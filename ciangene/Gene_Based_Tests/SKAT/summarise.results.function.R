suppressPackageStartupMessages(library(snpStats))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(xtable))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(rentrez))

ldak='/cluster/project8/vyp/cian/support/ldak/ldak'
source('/SAN/vyplab/UCLex/scripts/DNASeq_pipeline/ciangene/Gene_Based_Tests/plot/annotate_qqplot.R') ## Modified qqplot so I can label specific genes. 
lof <-  c("frameshift deletion", "frameshift substitution", "frameshift insertion",  "stoploss SNV", "splicing"
		,"stopgain SNV","exonic;splicing"
		)
snps<-read.table("/SAN/vyplab/UCLex/mainset_September2016/cian/Annotations/func.tab",sep='\t',header=T) 
snps.filt<-snps[snps$ExonicFunc %in% lof & snps$ExAC_MAF < 0.0001,]
snps.filt<-snps.filt[!duplicated(snps.filt$SNP),]

data<-'/SAN/vyplab/UCLex/mainset_September2016/cian/allChr_snpStats_out'
Rscript='/cluster/project8/vyp/vincent/Software/R-3.3.0/bin/Rscript'
gviz.script='/SAN/vyplab/UCLex/scripts/DNASeq_pipeline/ciangene/Gene_Based_Tests/SKAT/gviz.R'

bamList<-read.table('/SAN/vyplab/UCLex/support/SampleBamList',sep='\t')

summarise<-function(dir,genes=NULL,outputDirectory='Results',plot=TRUE,Title=basename(outputDirectory),percent=5,disease='heart',gviz=TRUE)
{
	dir<-paste0(dir,'/')
	outputDirectory<-paste0(outputDirectory,'/')
#	if(file.exists(outputDirectory)) file.remove(outputDirectory)
	if(!file.exists(outputDirectory)) dir.create(outputDirectory,recursive=TRUE)

	file.copy('/SAN/vyplab/UCLex/support/SKAT/README.docx',paste0(outputDirectory,'README.docx')) 
	system(paste('chmod 777',paste0(outputDirectory,'README.docx')))
	files<-list.files(dir,pattern='csv',full.names=T,recursive=T)
	files<-files[grep("chr",files)]
	dirName<-paste(basename(dirname(dirname(dir))),basename(dirname(dir)),basename(dir))
	if(length(files)!=23)message( paste('Problem-',dirName,'might be missing a few chromosomes...' ) )

	#file.copy(paste0( dirname(files)[1],'/qc/case_pca.plot.pdf'),paste0(outputDirectory,'Case_PCA_plots.pdf')) 
	case.qc.file<-paste0( dirname(files)[1],'/qc/cases_removed_because_of_low_read_depth.tab')
	if(file.size(case.qc.file)>0)file.copy(case.qc.file,paste0(outputDirectory,'cases_removed_because_of_low_read_depth.tab')) 
	#file.copy(paste0( dirname(files)[1],'/case.list'),paste0(outputDirectory,'case.list.txt')) 

	message("Added PCA plot of cases to outputDirectory")
	oFile<-paste0(outputDirectory,Title,'_SKAT_results.csv') 
	if(file.exists(oFile))file.remove(oFile)
	message("Merging results found in:")
	for(i in 1:length(files))
	{
		print(files[i])
		if(file.size(files[i])>0)
		{
			file<-read.csv(files[i],header=T)
			names<-colnames(file)
			write.table(file,oFile,col.names=!file.exists(oFile),row.names=F,quote=T,sep=',',append=T)
		}
	}

	file<-read.csv(oFile,header=F)
	file.remove(oFile)
	colnames(file)<-names
	file$SKATO<-as.numeric(as.character(file$SKATO))
	file<-file[order(file$SKATO),]
	file<-file[-grep("nb.cases",file$nb.cases),]
	#write.table(file,oFile,col.names=T,row.names=F,quote=T,sep=',')

	MergeFiles<-function(files,out)
	{
		if(file.exists(out))file.remove(out)
		for(carry in 1:length(files))
		{
			if(file.size(files[carry])>0)
			{
				c.file<-unique(read.csv(files[carry],header=F,sep='\t'))
				write.table(c.file, out,col.names=F,row.names=F,quote=T,sep=',',append=T)
			}
		}
		outFile<-read.csv(out,header=FALSE)
		return(outFile)
	}

	message('Merging case carriers')
	carriers<-list.files(dir,pattern='case_carriers',full.names=T,recursive=T)
	carry.oFile<-paste0(outputDirectory,Title,'_case_carriers.csv') 
	if(file.exists(carry.oFile))file.remove(carry.oFile)
	carriers<-MergeFiles(carriers,carry.oFile)
	message('Merging ctrl carriers')

	ctrl.carriers<-list.files(dir,pattern='ctrl_carriers',full.names=T,recursive=T)
	ctrl.carry.oFile<-paste0(outputDirectory,Title,'_control_carriers.csv') 
	file.remove(ctrl.carry.oFile)
	ctrl.carriers<-MergeFiles(ctrl.carriers,ctrl.carry.oFile)
	message('Merging results by SNP')
	if(file.exists(ctrl.carry.oFile))file.remove(ctrl.carry.oFile)

	by.snp<-list.files(dir,pattern='SKAT_results_by_SNP',full.names=T,recursive=T)
	by.snp.out<-paste0(outputDirectory,Title,'_results_by_SNP_cases.csv') 
	by.snp<-MergeFiles(by.snp,by.snp.out)
	colnames(by.snp)<-c("SNP", "case.snp.hets", "case.snp.homs", "case.mafs.snp",
		"ctrl.snp.hets", "ctrl.snp.homs", "ctrl.mafs.snp", "ENSEMBL",
		"Symbol", "SKATO", "nb.snps", "nb.cases", "nb.ctrls", "nb.alleles.cases",
		"nb.alleles.ctrls", "case.maf", "ctrl.maf", "total.maf", "nb.case.homs",
		"nb.case.hets", "nb.ctrl.homs", "nb.ctrl.hets", "Chr", "Start",
		"End", "FisherPvalue", "OddsRatio", "CompoundHetPvalue", "minCadd",
		"maxExac", "min.depth", "MeanCallRateCases", "MeanCallRateCtrls",
		"MaxMissRate", "HWEp", "MinSNPs", "MaxCtrlMAF", "SNPs", "GeneRD",
		"CaseSNPs", "SKATbeSNPs")
	message('Merging compound hets files cases')

	cpd.ctrls<-FALSE
	cpd.cases<-FALSE
	CompoundHets_cases<-list.files(dir,pattern='CompoundHets_cases',full.names=T,recursive=T)
	if(length(CompoundHets_cases)>0)
	{
	CompoundHets_cases.oFile<-paste0(outputDirectory,Title,'_CompoundHets_cases.csv') 
	CompoundHets_cases<-MergeFiles(CompoundHets_cases,CompoundHets_cases.oFile)
	message('Merging compound hets files ctrls')
	cpd.cases<-TRUE
	}

	#CompoundHets_ctrls<-list.files(dir,pattern='CompoundHets_ctrls',full.names=T,recursive=T)
	#if(length(CompoundHets_ctrls)>0)
	#{
	#CompoundHets_ctrls.oFile<-paste0(outputDirectory,Title,'_CompoundHets_ctrls.csv') 
	#CompoundHets_ctrls<-MergeFiles(CompoundHets_ctrls,CompoundHets_ctrls.oFile)
	#cpd.ctrls<-TRUE
	##if(file.exists(CompoundHets_ctrls.oFile))file.remove(CompoundHets_ctrls.oFile)
	#}

	SummariseCpdHets<-function(file)
	{
		file$LOF.rare.snps<-0 
		for(i in 1:nrow(file))
		{
			file$LOF.rare.snps[i]<-length(which(unlist(strsplit(file$V3[i],';')) %in% snps.filt$SNP))
		}
		filt<-subset(file,file$LOF.rare.snps>1)
		colnames(filt)<-c("Gene",'Sample','SNPs','Nb.minor.alleles.total','nb.lof.rare.snps')
		return(filt)
	}
	message('Processing cpd hets...')
	if(cpd.cases)
	{
		cpd.hets.cases<-SummariseCpdHets(CompoundHets_cases)
		if(nrow(cpd.hets.cases)>0)write.table(cpd.hets.cases,CompoundHets_cases.oFile,col.names=T,row.names=F,quote=F,sep='\t') else file.remove (CompoundHets_cases.oFile)
	}
	#if(cpd.ctrls)
	#{
	#	cpd.hets.ctrls<-SummariseCpdHets(CompoundHets_ctrls)
	#	if(nrow(cpd.hets.ctrls)>0)write.table(cpd.hets.ctrls,CompoundHets_ctrls.oFile,col.names=T,row.names=F,quote=F,sep='\t') else file.remove (CompoundHets_ctrls.oFile)
	#}

	message("Making list of samples that are carriers per variant")
	file$Carriers<-0
	for(snp in 1:nrow(carriers))
	{
		if(!is.na(carriers[snp,1]))
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
	}

	file$Nb.Carriers<-0 ## I dont want signal to be driven by a small number of cases which would prob be an artefact
	file$Nb.case.snps<-0 ## I dont want signal to be driven by a small number of cases which would prob be an artefact
	for(car in 1:nrow(file))
	{
		file$Nb.Carriers[car]<-length(unique(unlist(strsplit(as.character(file$Carriers[car]),';')) ))
		file$Nb.case.snps[car]<-length(unique(unlist(strsplit(as.character(file$CaseSNPs[car]),';')) ))
	}
	if(length(grep('ADA',colnames(file)))>0 )  file$ADA<-NULL
	write.table(file,paste0(outputDirectory,Title,'_SKAT.csv'),col.names=T,row.names=F,quote=T,sep=',',append=F)


############################################################################################################
############################################################################################################
	## Now make a filtered list of more plausible results
	if(!is.null(percent))
	{
		percent.cases.carriers<-percent
		nb.cases.required<- round( as.numeric(file$nb.cases) * (percent.cases.carriers/100)  ) [1]
	} else 
	if(mean(file$nb.cases)>20) nb.cases.required<- 5 else nb.cases.required <-2
	file<-file[which(file$SKATO>0),]
	if(!is.null(genes))
	{
		if( file.exists(file.path(genes))) genes<-read.table(genes,header=FALSE)[,1]
		file$Candidate<-FALSE
		file$Candidate[file$Symbol %in% genes[,1]]<-TRUE
		file<-data.frame(cbind(file[,ncol(file)],file[,1:(ncol(file)-1)]))
		colnames(file)[1]<-'Candidate'
	}# else{
	pval<-0.000001
	filt<-subset(file,file$Nb.Carriers>=nb.cases.required & file$MeanCallRateCases >0.8 & file$MeanCallRateCtrls > 0.8 & file$Nb.case.snps > 2) 
	if(nrow(filt)==0) filt<-subset(file,file$Nb.Carriers>=nb.cases.required & file$MeanCallRateCases >0.8 & file$MeanCallRateCtrls > 0.8) 
	
	plot.file<-filt 
	filt$FisherPvalue<-as.numeric(filt$FisherPvalue)
	filt<-subset(filt, as.numeric(filt$CompoundHetPvalue)<=pval | ( as.numeric(filt$SKATO)<=pval | filt$FisherPvalue<=pval))


	if(nrow(filt)>0)
	{

		# make a table of most likely causative gene for each case. 
		cases<-unique(unlist(strsplit(filt$Carriers,split=';') ))
		case.dat<-data.frame(matrix(nrow=length(cases),ncol=5))
		case.dat[,1]<-cases
		colnames(case.dat)<-c('Case','CausativeGene','Variants','RareAlleleCount','Phenopolis')
		for(case in 1:length(cases))
		{
			case.genes<-filt$Symbol[grep(cases[case],filt$Carriers)][1] # just take first gene
			for(gen in 1:length(case.genes))
			{
				gene.variants<-carriers[carriers$V3==filt$ENSEMBL[filt$Symbol==case.genes[gen]],]
				case.gene.variants<-paste(gene.variants[gene.variants$V2==cases[case],1],collapse=',')
				if(gen==1)case.variants<-case.gene.variants else case.variants<-c(case.variants,case.gene.variants)
				allele.count<-paste(gene.variants[gene.variants$V2==cases[case],5],collapse=',')
				if(gen==1)allele.count.all<-allele.count else allele.count.all<-c(allele.count,allele.count.all)
			}
			case.variants<-paste(unlist(case.variants),collapse=';')
			case.dat$CausativeGene[case]<-case.genes
			case.dat$Variants[case]<-case.variants
			case.dat$RareAlleleCount[case]<-allele.count.all
			phen<-paste0('https://uclex.cs.ucl.ac.uk/variant/',gsub(case.dat$Variants[1],pattern='_',replacement='-') )
			case.dat$Phenopolis[case]<-phen

		}
		write.table(case.dat,paste0(outputDirectory,Title,'_solved_cases.csv'),col.names=T,row.names=F,quote=T,sep=',',append=F)

		filt$Nb.relevant.papers<-0
		message('Finding relevant papers...')
		for(ro in 1:nrow(filt))
		{
			ff<- entrez_search(db="pubmed", term= paste('(',filt$Symbol[ro],')AND(',disease,')'))
			filt$Nb.relevant.papers[ro]<-ff$count
		}
		filt$pubmed.disease.term<-disease # record in output what disease we searched pubmed for with gene name
		SKATout<-paste0(outputDirectory,Title,'_SKAT_filtered.csv')
		write.table(filt,SKATout,col.names=T,row.names=F,quote=T,sep=',',append=F)

		if(nrow(filt)>10)filt<-filt[1:10,]
		rownames(filt)<-1:nrow(filt)

		message("Making HTML table for top genes")
		filt.xtable<-xtable(filt,caption=paste(Title,"SKAT top genes") ,digits=2, display = c(rep("s",4),'E',rep("d",5),rep("E",6),rep('d',3),rep('E',11),rep("s",7),rep('d',3),rep('s',1)))
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


		RData.file<-paste0(outputDirectory,Title,'_prep.RData')
		message("Saving workspace to ",RData.file)
		save(list=ls(environment()),file=RData.file)

		if(gviz)
		{
			message('Now plotting top genes/samples')
			variantDir<-file.path(outputDirectory,'topGenePlots')
		 	if(!file.exists(variantDir))dir.create(variantDir,recursive=TRUE)
			for(gene in 1:nrow(filt))
			{
				carriers<-unlist(strsplit(filt$Carriers[gene],';') )
				nb.carriers<-length(carriers)
				if(nb.carriers>20) nb.samples<-20 else nb.samples<-nb.carriers
				for(carrier in 1:nb.samples)
				{
					sample<-carriers[carrier]
					bams<-bamList[bamList[,2] %in% sample,]
					bams<-bams[grep('sorted',bams[,1]) ,]
					if(nrow(bams)>1)bams<-bams[which(file.size(bams[,1]) ==  max(file.size(bams[,1]) )),]

					pdf<-paste0(variantDir,'/',filt$Symbol[gene],'_',sample,'.pdf')
			#		run<-paste(Rscript,gviz,'--skat',SKATout,'--outPDF',pdf, '--Sample',sample,'--Gene', filt$Symbol[gene],'--sampleBam', bams[,1])
					run<-paste(Rscript,gviz.script,'--skat',SKATout,'--outPDF',pdf, '--Sample',sample,'--Gene', filt$Symbol[gene]) # bam reads display not clean yet
					message(run)
					system(run)
				}
			}
		}
	

	#}

	if(plot & (nrow(plot.file)>1000)  ) 
	{
		file<-plot.file
		qqplot.name<-paste0(paste0(outputDirectory,'qqplots.pdf')) 
		message(paste('QQplot stored in:',qqplot.name))
		pdf(qqplot.name)

		if(!is.null(genes)) 
		{			
			tst <- qq.chisq(-2*log(file$SKATO), df=2, pvals=TRUE,main=Title)
			labelQQplot(pval.labels=file$Symbol,qqplotOut=tst,pvals.raw=file$SKATO, genelist= genes) 
		} else{
			sig.genes<-subset(gene.pvals,gene.pvals$SKATO <0.000001)$Symbol
			if(length(sig.genes)>5)sig.genes<-sig.genes[1:5]
			tst <- qq.chisq(-2*log(file$SKATO), df=2,pvals=TRUE,main=Title)
			labelQQplot(pval.labels=file$Symbol,qqplotOut=tst,pvals.raw=file$SKATO, genelist= sig.genes )
		}
		dev.off()
	}#plot

	}#nrow(filt)>0)

	#message("Saving workspace to ",RData.file)
	#save(list=ls(environment()),file=RData.file)
} # summarise
