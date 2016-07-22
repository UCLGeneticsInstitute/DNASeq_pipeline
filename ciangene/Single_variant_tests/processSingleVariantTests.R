library(snpStats)

dir<-'/scratch2/vyp-scratch2/cian/UCLex_February2015/Single_variant_tests'
files<-list.files(dir,full.names=T)
tests<-c('counts_assoc.qassoc$','fisher.qassoc$','perm.assoc.fisher','tk_depth.assoc.linear$')

cohort.list<-c('Levine','Hardcastle','IoO','IoN','Kelsell','LambiaseSD','Lambiase_','LayalKC','Nejentsev','PrionUnit','Prionb2','Shamima','Sisodiya','Syrris','Vulliamy','WebsterURMD'
)

dat<-data.frame(matrix(nrow=length(cohort.list),ncol=(length(tests)+1)))
dat[,1]<-cohort.list
colnames(dat)<-c('Pheno',tests)

oFile<-'SingleVariantTests.pdf'
pdf(oFile)
par(mfrow=c(2,2))
for(i in 1:length(cohort.list))
{
	cohort<-cohort.list[i]
	cohort.files<-files[grep(cohort,files)]
	for(test in 1:length(tests))
	{	
		inFile<-cohort.files[grep(tests[test],cohort.files)]
		system(paste('tr -s " " < ',inFile,'>',paste0(inFile,'_clean')))
		print(inFile)
		file<-read.table(paste0(inFile,'_clean'),header=T)
		main=paste(cohort,tests[test])
		qq.chisq(-2*log(file$P),df=2,x.max=30,main=main)
	}
}

dev.off()