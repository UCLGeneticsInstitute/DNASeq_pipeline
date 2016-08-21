## changes
# added fisher test and odds ratio. 
# changed nb.cases/ctrls to nb non NA calls (so nb.cases now is nb.patients*nb.variants)
# changed case.count/ctlr.count to allele counts. 1 het and 1 hom is a count of 3 for nb.variants.cases 
# not removing non caucasians - now adding top two PCs as covariates into SKATO instead. 
# removing related individuals. 
# Fixed snp/gene matching
# Keeping only damaging variants

library(SKAT)
library(snpStats)
library(HardyWeinberg)
library(stringr)

oFile<-paste0(oDir,pheno.matching[phen,1],'_snps')
if(file.exists(oFile))file.remove(oFile)

caseFile<-paste0(oDir,pheno.matching[phen,1],'_case_carriers')
if(file.exists(caseFile))file.remove(caseFile)
ctrlFile<-paste0(oDir,pheno.matching[phen,1],'_ctrl_carriers')
if(file.exists(ctrlFile))file.remove(ctrlFile)

results$nb.variants.cases<-0
results$nb.variants.ctrls<-0
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
results$FisherPvalue<-0
results$OddsRatio<-0

#results$nb.unique.func.variants.cases<-0
#results$nb.unique.lof.variants.cases<-0
#results$nb.unique.func.variants.ctrls<-0
#results$nb.unique.lof.variants.ctrls<-0

### columns below this arent working yet. 
#results$nb.compound.func.hets.cases<-0
#results$nb.compound.func.hets.cases<-0
#results$nb.compound.func.hets.ctrls<-0
#results$nb.compound.lof.hets.cases<-0
#results$nb.compound.lof.hets.ctrls<-0


func<-read.csv(paste0(dirname(oDir),'/Annotations/func.tab'),header=TRUE,sep='\t')
func$isFunc<-FALSE
func$isFunc[func$Func&&func$Rare]<-TRUE
func$isLOF<-FALSE
func$isLOF[func$LOF&&func$Rare]<-TRUE

write.table( data.frame( colnames(clean.pheno.snps), current.pheno),paste0(oDir,pheno.matching[phen,1],'_case_control_list'), col.names=F,row.names=F,quote=F,sep='\t') 

#results<-results[with(results, order(Chr, Start,End)), ]
uniq.genes<-unique(results[,1])
nb.genes<-length(uniq.genes)

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

	results$Chr[gene]<-gene.chr
	results$Start[gene]<-gene.start
	results$End[gene]<-gene.end

	if(nb.snps.in.gene>0)
	{
		test.gene<-unique(snp.gene$ENSEMBL[snp.gene$SNP %in% rownames(gene.snp.data)])  ## match to uniq. genes as a check, 
		if(length(unique(test.gene))>1) message ("SNPs span multiple genes")
		if(test.gene!=uniq.genes[gene]) stop ("Genes not sorted correctly")
		if(results$ENSEMBL[gene]!=uniq.genes[gene]) stop ("Genes not sorted correctly")

		case.snps<-gene.snp.data[,current.pheno==1]
		ctrl.snps<-gene.snp.data[,current.pheno==0]

		maf.snp.cases<-apply(case.snps,1, function(x) signif(maf(as.numeric(unlist(table(unlist(x))))),2)  )
		maf.snp.ctrls<-apply(ctrl.snps,1, function(x) signif(maf(as.numeric(unlist(table(unlist(x))))),2)  )

		damaging.snps<-names(which(maf.snp.cases>maf.snp.ctrls))
		final.snp.set<-gene.snp.data[rownames(gene.snp.data) %in% damaging.snps, ]
		nb.snps.in.gene2<-nrow(final.snp.set)
		
		if(nb.snps.in.gene2>0)
		{
			case.snps<-final.snp.set[,current.pheno==1]
			ctrl.snps<-final.snp.set[,current.pheno==0]
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
			results$SKATO[gene] <- SKAT(t(as.matrix(final.snp.set)) , obj, missing_cutoff=0.8, estimate_MAF=2,method="optimal.adj")$p.value
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
				case.car<- rownames( data.frame(unlist(apply(case.snps,1,function(x) which(x>0 )))))
				if(nrow(case.snps)==1)
				{
					case.carriers<-colnames(case.snps)[apply(case.snps,1,function(x) which(x>0 ))]
					case.variants<-str_extract(rownames(case.snps),"[0-9]{1,2}_[0-9]+_[A-Z]_[A-Z]")
					case.dat<-data.frame(cbind(case.variants,case.carriers))

				} else
				{
					case.carriers<-case.car
					case.carriers.clean<-gsub(case.carriers,pattern="[0-9]{1,2}_[0-9]+_[A-Z]_[A-Z]\\.",replacement="")
					case.variants<-str_extract(case.car,"[0-9]{1,2}_[0-9]+_[A-Z]_[A-Z]")
					case.dat<-data.frame(cbind(case.variants,case.carriers.clean))
				}
				
				if( ( identical(case.dat[,1],case.dat[,2])  | is.na(case.dat[,1])) && nrow(case.snps)>1)
				{
					case.carriers<-colnames(case.snps)[apply(case.snps,1,function(x) which(x>0 ))]
					case.variants<-str_extract(rownames(case.snps),"[0-9]{1,2}_[0-9]+_[A-Z]_[A-Z]")
					case.dat<-data.frame(cbind(case.variants,case.carriers))
				} 
				if( ( identical(case.dat[,1],case.dat[,2])  | is.na(case.dat[,1]) ) && nrow(case.snps)==1)
				{
					case.carriers<-case.car
					case.carriers.clean<-gsub(case.carriers,pattern="[0-9]{1,2}_[0-9]+_[A-Z]_[A-Z]\\.",replacement="")
					case.variants<-str_extract(case.car,"[0-9]{1,2}_[0-9]+_[A-Z]_[A-Z]")
					case.dat<-data.frame(cbind(case.variants,case.carriers.clean))
				} 
			case.dat
			}

			case.dat<-fixNames(case.snps)
			ctrl.dat<-fixNames(ctrl.snps)
			case.dat<-data.frame(case.dat,uniq.genes[gene])
			ctrl.dat<-data.frame(ctrl.dat,uniq.genes[gene])
			write.table(case.dat, caseFile, col.names=FALSE,row.names=F,quote=F,sep='\t',append=T)
			write.table(ctrl.dat,ctrlFile, col.names=FALSE,row.names=F,quote=F,sep='\t',append=T)


			## Now do compound hets
			func.snps<-subset(func,isFunc)

			gene.func<-final.snp.set[ rownames(final.snp.set) %in% func.snps$SNP,]

			gene.lof<-final.snp.set[ rownames(final.snp.set) %in% func.snps$SNP,]



			#results$nb.compound.func.hets.cases
			#results$nb.compound.func.hets.ctrls
			#results$nb.compound.lof.hets.cases
			#results$nb.compound.lof.hets.ctrls

		}
	}
}

results <- results[order(results[,2]), ]
results.out <- paste0(oDir,pheno.matching[phen,1],'_skat.csv')
results<-results[order(results$SKATO),]
qqplot.out <- paste0(dirname(oDir),'/plots/',pheno.matching[phen,1],'_skat_QQplot.png')

write.table(results,results.out,col.names=T,row.names=F,quote=F,sep=',')
png(qqplot.out)
qq.chisq(-2*log(results$SKATO), df=2, x.max=30, main=paste(pheno.matching[phen,1],'SKAT'),cex.main=0.7)	
dev.off() 



