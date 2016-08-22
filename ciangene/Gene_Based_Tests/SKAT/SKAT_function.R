## how to use me. 
#  Rscript='/share/apps/R/bin/Rscript'
#  SKAT=/SAN/vyplab/UCLex/scripts/DNASeq_pipeline/ciangene/Gene_Based_Tests/SKAT/SKAT_function.R
#  CaseFile = path to one column file containing list of cases
#  ControlFile = path to one column file containing list of cases. If not specified, will use all non cases. 
#  $Rscript $SKAT --case.list $CaseFile --oDir outputDirectory --control.list $ControlFile

### Changes ####
# Added CADD and ExAC maf filter. 
# Added ability to specify minReadDepth for filtering - default zero is fastest. 
# Added ability to specify genes you're interested in testing: --TargetGenes for a file with names or --SampleGene for a single gene name
# Add CompoundHeterozygote Function and pvalue.  
# changed missingess to missing_cutoff=0.2. So variants with missingess >20% will be removed. 
# added fisher test and odds ratio. 
# changed nb.cases/ctrls to nb non NA calls (so nb.cases now is nb.patients*nb.variants)
# changed case.count/ctlr.count to allele counts. 1 het and 1 hom is a count of 3 for nb.variants.cases 
# now not removing non caucasians - adding top two PCs as covariates into SKATO instead. 
# removing related individuals. 
# Fixed snp/gene matching
# Keeping only damaging variants
########

ldak='/cluster/project8/vyp/cian/support/ldak/ldak'

suppressPackageStartupMessages(library("SKAT"))
suppressPackageStartupMessages(library("snpStats"))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("HardyWeinberg"))
suppressPackageStartupMessages(library("stringr"))

option_list <- list(
	make_option(c("-v", "--verbose"), action="store_true", default=TRUE,help="Print extra output [default]"),
 	make_option(c("--case.list"),  help="one column file containing list of cases",type='character'),
 	make_option(c("--control.list"), default=NULL, help="one column file containing list of controls",type='character'),
 	make_option(c("--SampleGene"), default=NULL, help="Gene Symbol",type='character'),
 	make_option(c("--TargetGenes"), default=NULL, help="Gene File",type='character'),
 	make_option(c("--MinReadDepth"), default=0, help="Specify MinReadDepth",type='character'),
 	make_option(c("--SavePrep"), default=FALSE, help="Do you want to save an image of setup?",type='character'),
 	make_option(c("--minCadd"), default=20, help="minimum CADD score for retained variants",type='character'),
 	make_option(c("--maxExac"), default=0.01, help="Max EXAC maf for retained variants",type='character'),
 	make_option(c("--oDir"),default='SKATtest',type='character')
 )


opt <- parse_args(OptionParser(option_list=option_list))
if ( opt$verbose ) {
 write("Starting Argument checks...\n", stderr())
}


case.list=opt$case.list
control.list=opt$control.list
outputDirectory=opt$oDir

if( (!is.null(opt$SampleGene)) && ( !is.null(opt$TargetGenes)) ) stop("Please supply either a single gene to SampleGene or a file of gene names to TargetGenes. Not both.")
if(!is.null(opt$SampleGene))
{
 	TargetGenes=opt$SampleGene
 	message(paste("Setting up environment to test for gene:",TargetGenes))
} 
if(!is.null(opt$TargetGenes))
{
	TargetGenes=read.table(opt$TargetGenes,header=FALSE)[,1]
	message(paste("Read",length(TargetGenes),'genes from file',opt$TargetGenes))
} 

min.depth=as.numeric(opt$MinReadDepth) 
SavePrep=opt$SavePrep
maxExac=as.numeric(opt$maxExac)
minCadd=as.numeric(opt$minCadd)

## Check cases
if(!file.exists(case.list))stop("Case file doesn't exist. ")
message(paste('Reading cases from',paste0('--',case.list,'--')))
case.list<-read.table(case.list,header=FALSE)
message(paste('Found',nrow(case.list),'cases'))
case.list<-case.list[,1]


## Set up controls
if(is.null(control.list))message("Controls not specified, so will use all non case samples as controls.")
if(!is.null(control.list))
{
	if(!file.exists(control.list))message("Control file doesn't exist, so will use all non case samples as controls.")
	message(paste('Reading controls from',control.list))
	control.list<-read.table(control.list,header=FALSE)
	message(paste('Found',nrow(control.list),'cases'))
	control.list<-control.list[,1]
}

## make or break oDir
if(outputDirectory=='SKATtest')message('output directory not specified so will use default folder of ./SKATtest')


message("Finished Argument checks.\n")
message("Starting test setup.\n")

doSKAT<-function(case.list,control.list=NULL,outputDirectory,min.depth=0,release='July2016',compoundHets='Yes',TargetGenes=NULL,
	minCadd=20,maxExac=0.01
	)
{
	outputDirectory<-paste0(outputDirectory,'/')
	rootODir<-paste0('/SAN/vyplab/UCLex/mainset_',release,'/cian/') 
	if(!file.exists(outputDirectory)) dir.create(outputDirectory)
	qc<-paste0(outputDirectory,'/qc/')
	if(!file.exists(qc))dir.create(qc)

	oFile<-paste0(outputDirectory,'SKAT_results_by_SNP.tab')
	if(file.exists(oFile)) file.remove(oFile)
	caseFile<-paste0(outputDirectory,'case_carriers')
	if(file.exists(caseFile)) file.remove(caseFile)
	ctrlFile<-paste0(outputDirectory,'ctrl_carriers')
	if(file.exists(ctrlFile)) file.remove(ctrlFile)

	message("Reading in snp/gene database")
	snp.gene.base<-read.table(paste0(rootODir,'snp_gene_id'),header=F)
	colnames(snp.gene.base)<-c("SNP",'ENSEMBL')
	gene.dict<-read.table(paste0(rootODir,'gene_dict_skat'),header=T)
	gene.dict$ENST<-NULL
	gene.dict<-unique(gene.dict)
	snp.gene<-merge(gene.dict,snp.gene.base,by="ENSEMBL",all.y=T)
	
	## This anno file contains Exac mafs, CADD scores etc. Filter SNP list based on user input. 
	snp.annotations<-read.table(paste0(rootODir,'Annotations/func.tab'),header=TRUE,sep='\t')
	filtered.snp.list<-subset(snp.annotations,snp.annotations$CADD>=minCadd & snp.annotations$ExAC_MAF<=maxExac)$SNP
	message(paste(length(filtered.snp.list),'SNPs kept after CADD and Exac filters of',minCadd,'and',maxExac,'respectively.')) 

	if(!is.null(TargetGenes)) ## if we're testing a couple genes only, Im subsetting SNPmatrix to make it faster to read in. 
	{
		data<-paste0(rootODir,'gene_test_variants_out') 
		target.snp.info<- snp.gene[ snp.gene$Symbol %in% TargetGenes,]
		target.snps<- target.snp.info$SNP[target.snp.info$SNP %in% filtered.snp.list]
		write.table(target.snps,paste0(outputDirectory,TargetGenes[1]),col.names=F,row.names=F,quote=F,sep='\t')
		message(paste(length(target.snps),'SNPs found in:', TargetGenes))
		message("Extracting from dataset...")
		run<-paste(ldak,'--make-sp',paste0(outputDirectory,TargetGenes),'--sp',data,'--extract',paste0(outputDirectory,TargetGenes[1]))
		print(run)
		system(run)
		message("Reading in SNP data")
		snp.data<-paste0(outputDirectory,TargetGenes[1],'_out')
		snp.sp<-read.table(paste0(outputDirectory,TargetGenes[1],'_out.sp'),header=F)
		snp.bim<-read.table(paste0(outputDirectory,TargetGenes[1],'_out.bim'),header=F)
		snp.fam<-read.table(paste0(outputDirectory,TargetGenes[1],'_out.fam'),header=F)
		colnames(snp.sp)<-snp.fam[,1]
		rownames(snp.sp)<-snp.bim[,2]
	}

	if(is.null(TargetGenes))
	{
		write.table(filtered.snp.list,paste0(outputDirectory,'SNPlist'),col.names=F,row.names=F,quote=F,sep='\t')
		message("Extracting from dataset...")
		run<-paste(ldak,'--make-sp',paste0(outputDirectory,TargetGenes),'--sp',data,'--extract',paste0(outputDirectory,'SNPlist'))
		print(run)
		system(run)
		snp.sp<-read.table(paste0(rootODir,'gene_test_variants_out.sp'),header=F)
		snp.bim<-read.table(paste0(rootODir,'gene_test_variants_out.bim'),header=F)
		snp.fam<-read.table(paste0(rootODir,'gene_test_variants_out.fam'),header=F)
		colnames(snp.sp)<-snp.fam[,1]
		rownames(snp.sp)<-snp.bim[,2]
	}

	## ancestry matching
	ancestry<-read.table(paste0(rootODir,'UCLex_samples_ancestry'),header=T)
	#caucasians<-ancestry$V1[ancestry$Caucasian]
	#snps<-snp.sp[,colnames(snp.sp)%in%caucasians]
	## these are teh snps that have been filterd by func and maf
	unrelated<-read.table(paste0(dirname(rootODir),"/kinship/UCL-exome_unrelated.txt")) ## Keeping only unrelated individiuals
	snp.sp<-snp.sp[,colnames(snp.sp)%in%unrelated[,1]]
	print(paste('Read Depth filter set at:',min.depth ) ) 
	if(min.depth>0)
	{
		message("Reading in Read Depth data")
		read.depth<-read.table(paste0(rootODir,'Average.read.depth.by.gene.by.sample.tab'),header=T)
		sample.means<-colMeans(read.depth[,5:ncol(read.depth)])
		good.samples<-which(sample.means>=min.depth)
		bad.samples<-names(which(sample.means<min.depth))
		write.table(bad.samples,paste0(qc,'samples_removed_because_of_low_read_depth.tab'),col.names=F,row.names=F,quote=F,sep='\t')

		gene.means<-rowMeans(read.depth[,5:ncol(read.depth)])
		good.genes<-which(gene.means>=min.depth)
		bad.genes<-read.depth$Gene[which(gene.means<min.depth)]
		write.table(bad.genes,paste0(qc,'genes_removed_because_of_low_read_depth.tab'),col.names=F,row.names=F,quote=F,sep='\t')

		good.genes.snps<-snp.gene$SNP[snp.gene$ENSEMBL %in% read.depth$Gene[good.genes]]
		clean.snp.data<- snp.sp[rownames(snp.sp) %in% good.genes.snps ,  colnames(snp.sp) %in% colnames(read.depth)[good.samples] ]
		snp.gene.clean<-snp.gene[snp.gene$SNP %in%rownames(clean.snp.data),]
		write.table(snp.gene.clean,paste0(qc,'read.depth.ancestry.func.clean.snps.tab'),col.names=F,row.names=F,quote=F,sep='\t')

		good.genes.data<-snp.gene[snp.gene$ENSEMBL %in% read.depth$Gene[good.genes],]
		uniq.genes<-unique(good.genes.data$ENSEMBL)
		nb.genes<-length(uniq.genes)
		rm(read.depth)
	} else
	{
		message("No read depth filter specified. Skipping step.")
		clean.snp.data<- snp.sp
		good.genes.data<-snp.gene
		uniq.genes<-unique(good.genes.data$ENSEMBL)
		nb.genes<-length(uniq.genes)
	}

	###### clean
	rm(snp.sp,snp.gene.base)
	######

	current.pheno<-rep(NA,ncol(clean.snp.data))
	current.pheno[colnames(clean.snp.data)%in%case.list]<-1
	my.cases<-colnames(clean.snp.data)[colnames(clean.snp.data)%in%case.list]
	my.controls<-colnames(clean.snp.data)[colnames(clean.snp.data)%in%control.list]
	if(!is.null(control.list))
	{
		current.pheno[colnames(clean.snp.data)%in%control.list]<-0
	} else current.pheno[!colnames(clean.snp.data)%in%case.list]<-0
	write.table(my.cases,paste0(outputDirectory,'control.list'),col.names=FALSE,row.names=FALSE,quote=FALSE)
	write.table(my.controls,paste0(outputDirectory,'case.list'),col.names=FALSE,row.names=FALSE,quote=FALSE)

	## Make the outptu dataframe
	cols<-c("Gene",'SKATO','nb.snps','nb.cases','nb.ctrls','nb.variants.cases','nb.variants.ctrls','case.maf','ctrl.maf','total.maf','nb.case.homs',
		'nb.case.hets','nb.ctrl.homs','nb.ctrl.hets','Chr','Start','End','FisherPvalue','OddsRatio','CompoundHetPvalue','minCadd','maxExac'
		)
	results<-data.frame(matrix(nrow=nb.genes,ncol=length(cols)))
	colnames(results)<-cols
	results$Gene<-uniq.genes
	results$minCadd<-minCadd
	results$maxExac<-maxExacc
	results<-merge(gene.dict,results,by.y='Gene',by.x='ENSEMBL',all.y=T)
	srt<-data.frame(1:length(uniq.genes),uniq.genes)
	results<-merge(results,srt,by.y='uniq.genes',by.x='ENSEMBL')
	results<-results[order(results[,(ncol(results))]),]
	#results<-results[,1:(ncol(results)-1)]
	uniq.genes<-unique(results[,1])
	nb.genes<-length(uniq.genes)

	if(!is.null(TargetGenes))
	{
		results<-results[,1:(ncol(results)-1)]
		results<-results[results$Symbol %in% TargetGenes,]
		if(nrow(results)==0)stop('Specified genes not found ffs.')
		uniq.genes<-unique(results[,1])
		nb.genes<-length(uniq.genes)
	}

	if(SavePrep)
	{
		robj<-paste0(outputDirectory,'test_setup.RData')
		message(paste('Saving workspace image to', robj))
		save(list=ls(environment()),file=robj)
	}	
	###################################
	message(paste("Starting tests on", nb.genes, 'genes')) 
	###################################

	for(gene in 1:nb.genes)
	{
		if(gene%%100==0) message(paste('Up to gene',gene)) 
		gene.snps<-good.genes.data$SNP[ grep(uniq.genes[gene],good.genes.data$ENSEMBL) ]

		gene.data<- data.frame(t(data.frame(strsplit(gene.snps,'_')))) 
		gene.chr<-as.numeric(unique(gene.data[,1]))
		gene.start<-min(as.numeric(unique(gene.data[,2])))
		gene.end<-max(as.numeric(unique(gene.data[,2])))

		gene.snp.data<-clean.snp.data[ rownames(clean.snp.data) %in% gene.snps ,]
		nb.snps.in.gene<-nrow(gene.snp.data)
		print(paste(nb.snps.in.gene,'snps in', uniq.genes[gene],'before full filtering complete')) 

		results$Chr[gene]<-gene.chr
		results$Start[gene]<-gene.start
		results$End[gene]<-gene.end

		if(nb.snps.in.gene>0)
		{
			test.gene<-unique(snp.gene$ENSEMBL[snp.gene$SNP %in% rownames(gene.snp.data)])  ## match to uniq. genes as a check, 
			if(length(unique(test.gene))>1) message ("SNPs span multiple genes")
			if(test.gene!=uniq.genes[gene]) stop ("Genes not sorted correctly")
			if(results$ENSEMBL[gene]!=uniq.genes[gene]) stop ("Genes not sorted correctly")


			case.snps<-gene.snp.data[,which(current.pheno==1)]
			ctrl.snps<-gene.snp.data[,which(current.pheno==0)]

			maf.snp.cases<-apply(case.snps,1, function(x) signif(maf(as.numeric(unlist(table(unlist(x))))),2)  )
			maf.snp.ctrls<-apply(ctrl.snps,1, function(x) signif(maf(as.numeric(unlist(table(unlist(x))))),2)  )

			damaging.snps<-names(which(maf.snp.cases>maf.snp.ctrls))
			final.snp.set<-gene.snp.data[rownames(gene.snp.data) %in% damaging.snps, ]
			nb.snps.in.gene2<-nrow(final.snp.set)
			
			print(paste(nb.snps.in.gene2,'snps in', uniq.genes[gene],'after filtering'))
			if(nb.snps.in.gene2>0)
			{
				case.snps<-final.snp.set[,which(current.pheno==1)]
				ctrl.snps<-final.snp.set[,which(current.pheno==0)]
				results$nb.cases[gene]<-length(which(!is.na(unlist(case.snps))) )
				results$nb.ctrls[gene]<-length(which(!is.na(unlist(ctrl.snps))) )

				results$nb.variants.cases[gene]<-(length(grep(1,unlist(case.snps))))+ (length(grep(2,unlist(case.snps)))*2)
				results$nb.variants.ctrls[gene]<-(length(grep(1,unlist(ctrl.snps))))+ (length(grep(2,unlist(ctrl.snps)))*2)

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

				if(sum(final.snp.set,na.rm=T)>0)results$total.maf[gene]<-signif(maf(as.numeric(unlist(table(unlist(final.snp.set))))),2) 

				## these counts are for each snp in gene separately
				case.snp.hets<-apply(case.snps,1,function(x) length(grep(1,x)))
				case.snp.homs<-apply(case.snps,1,function(x) length(grep(2,x)))
				maf.snp.cases<-apply(case.snps,1, function(x) signif(maf(as.numeric(unlist(table(unlist(x))))),2)  )

				ctrl.snp.hets<-apply(ctrl.snps,1,function(x) length(grep(1,x)))
				ctrl.snp.homs<-apply(ctrl.snps,1,function(x) length(grep(2,x)))
				maf.snp.ctrls<-apply(ctrl.snps,1, function(x) signif(maf(as.numeric(unlist(table(unlist(x))))),2)  )

				ancestry.pcs<-ancestry[match(colnames(final.snp.set),ancestry$V1),]
				obj<-SKAT_Null_Model(current.pheno ~ ancestry.pcs$V3+ancestry.pcs$V4, out_type="D")
				#results$SKATO[gene] <- SKAT(t(as.matrix(gene.snp.data)) , obj, missing_cutoff=0.4, estimate_MAF=2)$p.value 
				results$SKATO[gene] <- SKAT(t(as.matrix(final.snp.set)) , obj, missing_cutoff=0.2, estimate_MAF=2,method="optimal.adj")$p.value
				results$nb.snps[gene] <- nrow(final.snp.set)


	       		mat<-matrix(c(results$nb.ctrls[gene]*2 - results$nb.variants.ctrls[gene],
	       						results$nb.variants.ctrls[gene],
	       						results$nb.cases[gene]*2 - results$nb.variants.cases[gene],
	       						results$nb.variants.cases[gene])
	       						, nrow = 2, ncol = 2)
				if (length(which(is.na(unlist(mat))))==0)
	       		{
	       			testy<-fisher.test(mat)
	       			results$FisherPvalue[gene]<-signif(testy$p.value,4) 
	       			results$OddsRatio[gene]<-signif(testy$estimate,4) 
	       		}


				snp.out<-data.frame( rownames(final.snp.set),case.snp.hets,case.snp.homs, maf.snp.cases,ctrl.snp.hets,ctrl.snp.homs, maf.snp.ctrls,results[gene,]) 
				write.table( snp.out, oFile, col.names=!file.exists(oFile),row.names=F,quote=F,sep='\t',append=T)

				fixNames<-function(snps)
				{
					car<- rownames( data.frame(unlist(apply(snps,1,function(x) which(x>0 )))))
					if(nrow(snps)==1)
					{
						carriers<-colnames(snps)[apply(snps,1,function(x) which(x>0 ))]
						variants<-str_extract(rownames(snps),"[0-9]{1,2}_[0-9]+_[A-Z]_[A-Z]")
						dat<-data.frame(cbind(variants,carriers))

					} else
					{
						carriers<-car
						carriers.clean<-gsub(carriers,pattern="[0-9]{1,2}_[0-9]+_[A-Z]_[A-Z]\\.",replacement="")
						variants<-str_extract(car,"[0-9]{1,2}_[0-9]+_[A-Z]_[A-Z]")
						dat<-data.frame(cbind(variants,carriers.clean))
					}
					
					if( ( identical(dat[,1],dat[,2])  | is.na(dat[,1])) && nrow(snps)>1)
					{
						carriers<-colnames(snps)[apply(snps,1,function(x) which(x>0 ))]
						variants<-str_extract(rownames(snps),"[0-9]{1,2}_[0-9]+_[A-Z]_[A-Z]")
						dat<-data.frame(cbind(variants,carriers))
					} 
					if( ( identical(dat[,1],dat[,2])  | is.na(dat[,1]) ) && nrow(snps)==1)
					{
						carriers<-car
						carriers.clean<-gsub(carriers,pattern="[0-9]{1,2}_[0-9]+_[A-Z]_[A-Z]\\.",replacement="")
						variants<-str_extract(car,"[0-9]{1,2}_[0-9]+_[A-Z]_[A-Z]")
						dat<-data.frame(cbind(variants,carriers.clean))
					} 
				return(dat)
				}


				case.dat<-fixNames(case.snps)
				ctrl.dat<-fixNames(ctrl.snps)
				case.dat<-data.frame(case.dat,uniq.genes[gene])
				ctrl.dat<-data.frame(ctrl.dat,uniq.genes[gene])
				write.table(case.dat, caseFile, col.names=FALSE,row.names=F,quote=F,sep='\t',append=T)
				write.table(ctrl.dat,ctrlFile, col.names=FALSE,row.names=F,quote=F,sep='\t',append=T)

				if(compoundHets=='Yes')
				{

					GetCompoundHets<-function(snp.dat,outFile)
					{
						compound.hets.names<-names(which(table(snp.dat$carriers.clean)>1))
						nb.hets<-0
						for(i in 1:length(compound.hets.names))
						{
							compound.snps<-t(snp.dat$variants[snp.dat$carriers.clean%in%compound.hets.names[i] ] )
							tt<-data.frame(compound.hets.names[i],compound.snps,uniq.genes[gene]) 
							if(length(tt)>0)
							{
								write.table(tt,outFile,col.names=F,row.names=F,quote=F,sep='\t',append=T)
								nb.hets<-nb.hets+1
							}
						}
						return(nb.hets)
					}
					compoundFileCases<-paste0(outputDirectory,'CompoundHets_cases')
					compoundFileCtrls<-paste0(outputDirectory,'CompoundHets_ctrls')

					caseTest<-FALSE
					ctrlTest<-FALSE
					if(length(unique(case.dat[,grep("carriers",colnames(case.dat))]))<nrow(case.dat) )
					{
						case.compound.hets<-GetCompoundHets(case.dat,compoundFileCases)
						caseTest<-TRUE
					}
					if(length(unique(ctrl.dat[,grep("carriers",colnames(ctrl.dat))]))<nrow(ctrl.dat) )
					{
						ctrl.compound.hets<-GetCompoundHets(ctrl.dat,compoundFileCtrls)
						ctrlTest<-TRUE
					}
					if(caseTest && ctrlTest)
					{
						nb.clean.cases<-ncol(case.snps)-case.compound.hets
						nb.clean.ctrls<-ncol(ctrl.snps)-ctrl.compound.hets

		       			mat<-matrix(c(nb.clean.ctrls,
		       						ctrl.compound.hets,
		       						nb.clean.cases,
		       						case.compound.hets)
		       						, nrow = 2, ncol = 2)
						results$CompoundHetPvalue[gene]<-fisher.test(mat,alternative='greater')$p.value
					}
				}


				print(results[gene,])
			}
		}
	}
###################################
 		results <- results[order(results[,2]), ]
		results.out <- paste0(outputDirectory,'skat.csv')
		results<-results[order(results$SKATO),]
		qqplot.out <- paste0(outputDirectory,'skat_QQplot.png')

		write.table(results,results.out,col.names=T,row.names=F,quote=F,sep=',')
	
	if(is.null(TargetGenes))
	{
		png(qqplot.out)
		qq.chisq(-2*log(results$SKATO), df=2, x.max=30, main='SKAT',cex.main=0.7)	
		dev.off() 
	}

print("Finished testing.")
} #doSKAT




#### now run function.
doSKAT(case.list=case.list,control.list=control.list,outputDirectory=outputDirectory,TargetGenes=TargetGenes,min.depth=min.depth,minCadd=minCadd,maxExac=maxExac)


