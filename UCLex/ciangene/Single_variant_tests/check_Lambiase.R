library(biomaRt)
bDir<-"/scratch2/vyp-scratch2/cian/UCLex_July2015/"
ldak<-"/cluster/project8/vyp/cian/support/ldak/ldak"
#load("/scratch2/vyp-scratch2/cian/UCLex_February2015/FastLMM_Single_Variant_all_phenos/dat.RData")
#write.table(dat$SNP,"Lambiase_snps",col.names=F,row.names=F,quote=F,sep="\t") 
source("/cluster/project8/vyp/cian/scripts/r/Filter_and_fisher.R")

file<-read.csv("/scratch2/vyp-scratch2/cian/UCLex_July2015/FastLMM_Single_Variant_all_phenos/Lambiase_final",header=T,sep="\t",stringsAsFactors=F)
file$Pvalue<-as.numeric(as.character(file$Pvalue))
file$TechKinPvalue<-as.numeric(as.character(file$TechKinPvalue))
dat<-subset(file,file$Pvalue<0.000001&file$TechKinPvalue>0.00001&file$C_A>3)
write.table(dat$SNP,"Lambiase_snps",col.names=F,row.names=F,quote=F,sep="\t") 
system(paste(ldak, "--make-sp snps.uclex --extract Lambiase_snps --sp", paste0(bDir, "allChr_snpStats")) )  

calls<-read.table("snps.uclex_out.sp",header=F) 
fam<-read.table("snps.uclex_out.fam",header=F) 
bim<-read.table("snps.uclex_out.bim",header=F) 
colnames(calls)<-fam[,1]
rownames(calls)<-bim[,2]

pheno<-read.table(paste0(bDir,'Clean_pheno_subset'))
fam<-read.table(paste0(bDir,'allChr_snpStats_out.fam'))
cohorts<-read.table(paste0(bDir,'cohort.list'))
colnames(pheno)<-c(rep("Samples",2),cohorts[,1])

calls<-calls[,colnames(calls)%in%pheno[,1] ] 
calls<-calls[,!is.na(pheno$Lambiase)]

cases<-grep("ambiase",colnames(calls))
case.group1<-grep("LambiaseSD",colnames(calls))
case.group2<-grep("Lambiase_April2014",colnames(calls))
case.group3<-cases[!cases %in% c(case.group1,case.group2) ]

countsCols<-c("SNP","CaseCounts","Case1.counts","Case2.counts","Case3.counts",
				"Nb.case1.nas","Nb.case2.nas","Nb.case3.nas","Ctrl.counts",
				"nb.ctrl.nas","nb.ctrls", 'nb.cases',
				"nb.case2.samples","nb.case1.samples","nb.case3.samples", "fisherP",'OR',"NoKinPvalue","TechKinPvalue",'ExtCtrl_MAF'
				) 
counts<-data.frame(matrix(nrow=nrow(dat),ncol=length(countsCols))) 
colnames(counts)<-countsCols
counts$SNP<-dat$SNP
counts$NoKinPvalue<-dat$Pvalue
counts$TechKinPvalue<-dat$TechKinPvalue
counts$ExtCtrl_MAF<-dat$ExtCtrl_MAF
for(i in 1:nrow(dat))
#for(i in 1:10)
{
	hit<-which(rownames(calls)%in%dat$SNP[i]) 
	if(length(hit)>1)stop(print("found at",i)) 
	case.counts<-sum(calls[rownames(calls)%in%dat$SNP[i],cases],na.rm=T) 	
	case.group1.counts<-sum(calls[rownames(calls)%in%dat$SNP[i],case.group1],na.rm=T) 	
	case.group2.counts<-sum(calls[rownames(calls)%in%dat$SNP[i],case.group2],na.rm=T) 	
	case.group3.counts<-sum(calls[rownames(calls)%in%dat$SNP[i],case.group3],na.rm=T) 	
	nb.case1.nas<-length(which(is.na(calls[rownames(calls)%in%dat$SNP[i],case.group1])) ) 
	nb.case2.nas<-length(which(is.na(calls[rownames(calls)%in%dat$SNP[i],case.group2])) ) 
	nb.case3.nas<-length(which(is.na(calls[rownames(calls)%in%dat$SNP[i],case.group3])) ) 

	ctrl.counts<-sum(calls[rownames(calls)%in%dat$SNP[i],-cases],na.rm=T) 	
	nb.nas.ctrls<-length(which(is.na(calls[rownames(calls)%in%dat$SNP[i],-cases])) )
	nb.ctrls<-length(which(!is.na(calls[rownames(calls)%in%dat$SNP[i],-cases])) )
	
	counts$CaseCounts[i]<-case.counts
	counts$Case1.counts[i]<-case.group1.counts
	counts$Case2.counts[i]<-case.group2.counts
	counts$Case3.counts[i]<-case.group3.counts
	counts$Nb.case1.nas[i]<-nb.case1.nas
	counts$Nb.case2.nas[i]<-nb.case2.nas
	counts$Nb.case3.nas[i]<-nb.case3.nas
	counts$nb.case1.samples[i]<-length(which(!is.na(calls[rownames(calls)%in%dat$SNP[i],case.group1])))
	counts$nb.case2.samples[i]<-length(which(!is.na(calls[rownames(calls)%in%dat$SNP[i],case.group2])))
	counts$nb.case3.samples[i]<-length(which(!is.na(calls[rownames(calls)%in%dat$SNP[i],case.group3])))
	counts$nb.cases[i]<-counts$nb.case1.samples[i]+counts$nb.case2.samples[i]+counts$nb.case3.samples[i]
	counts$Ctrl.counts[i]<-ctrl.counts
	counts$nb.ctrl.nas[i]<-nb.nas.ctrls
	counts$nb.ctrls[i]<-nb.ctrls
	
	mat<-matrix(c(counts$CaseCounts[i],
					(counts$nb.case1.samples[i]+counts$nb.case2.samples[i]+counts$nb.case3.samples[i])*2-counts$CaseCounts[i] , 
					counts$Ctrl.counts[i], 
					counts$nb.ctrls[i]*2-counts$Ctrl.counts[i] ),
		                       	nrow = 2, ncol = 2)
	if (length(which(is.na(unlist(mat))))==0)
	{
		counts$fisherP[i]<-fisher.test(mat)$p.value
		counts$OR[i]<-fisher.test(mat)$estimate 
	}

print(counts[i,])
}

save(counts,file="counts.RData") 

#fakes<- counts[counts$CaseCounts>0 & (counts$CaseCounts==counts$Case1.counts | counts$CaseCounts==counts$Case2.counts | counts$CaseCounts==counts$Case3.counts),]
#fakes<-fakes[order(-fakes$CaseCounts),]
#exit

keep<-subset(counts,counts$Case1.counts==CaseCounts|counts$Case2.counts==CaseCounts|counts$Case3.counts==CaseCounts)
keep<-subset(keep,keep$CaseCounts>0)
keep<-merge(keep,dat,by='SNP')

ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
filter="ensembl_gene_id"
attributes =  c("ensembl_gene_id", "external_gene_name")
#keep.anno <- getBM(attributes= attributes , filters = filter , values = keep$Gene , mart = ensembl)
#dat<-merge(keep,keep.anno,by.x='Gene',by.y='ensembl_gene_id',all.x=T)

genes<-read.table("/cluster/project8/vyp/cian/data/Cardio/Sudden_death/sudden_death/candidate_genes_neat_manual",header=F,sep="\t")
keep$Candidate<-FALSE
keep$Candidate[keep$external_gene_name%in%unlist(genes)]<-TRUE

calls<-prepData(keep)
pvals<-doFisher(calls,cases='Lambiase')
keep<-merge(keep,pvals,by='SNP')
write.table(keep,'keep',col.names=T,row.names=F,quote=F,sep="\t")


cases<-calls[rownames(cases)%in%file$SNP,grep("Lambiase",colnames(calls))]
for(i in 1:nrow(cases))
{
	dats<-t(data.frame(c(rownames(cases)[i],colnames(cases)[which(cases[i,]>0)]) ) ) 
	write.table(dats,'carriers',col.names=F,row.names=F,quote=F,sep="\t",append=T)
}


keep.anno <- getBM(attributes= attributes , filters = filter , values = rat$Gene , mart = ensembl)
rat<-merge(rat,keep.anno,by.x='Gene',by.y='ensembl_gene_id',all.x=T)
rat$Candidate<-FALSE
rat$Candidate[rat$external_gene_name%in%unlist(genes)]<-TRUE
rat$baseP<-as.numeric(as.character(rat$baseP))
rat$tkP<-as.numeric(as.character(rat$tkP))
rat$permeP<-as.numeric(as.character(rat$permeP))
rat<-rat[order(rat$tkP),]

by(rat$baseP,rat$Candidate,summary)
by(rat$tkP,rat$Candidate,summary)

tra<-subset(rat,rat$FILTER=='PASS')
tra<-subset(tra,tra$Call.rate>.8)
#tra<-subset(tra,tra$tkP<.0000001)


by(tra$baseP,tra$Candidate,summary)
by(tra$tkP,tra$Candidate,summary)


save(rat,dat,tra,file='Lambiase_single_variant_finding_false_positives.RData')



