library(GenABEL) 
library(ggplot2)
library(reshape) 

bDir<-paste0("/scratch2/vyp-scratch2/cian/UCLex_",release,'/')
dir<-paste0(bDir,"FastLMM_Single_Variant_all_phenos/") 
files<-list.files(dir,pattern="final",full.names=T) 
names<-gsub(basename(files),pattern="_.*",replacement="") 

getGIs<-function(file,pheno)
{

	file$Pvalue[file$Pvalue==0]<-0.0000001
	file$TechKinPvalue[file$TechKinPvalue==0]<-0.0000001
	file$PermP[file$PermP==0]<-0.0000001

	mafs<-c(0,0.00001,0.0001,0.001,0.01,0.1) 
	GIs<-data.frame(matrix(nrow=length(mafs),ncol=5) ) 
	colnames(GIs)<-c("Group","MAF","Base","TKRD","Perm") 
	GIs$MAF<-mafs
	GIs$Group<-pheno
	for(maf in 1:length(mafs))
	{
		dat<-subset(file,file$ExtCtrl_MAF>=mafs[maf]) 
	
		base<-dat$Pvalue[!is.na(dat$Pvalue)]
		tkrd<-dat$TechKinPvalue[!is.na(dat$TechKinPvalue)]
		perm<-dat$PermP[!is.na(dat$PermP)]

		GIs$Base[maf]<-estlambda(base)$estimate
		GIs$TKRD[maf]<-estlambda(tkrd)$estimate
		GIs$Perm[maf]<-estlambda(perm)$estimate
	}
GIs
}

for( i in 1:length(names))
{
	file<-read.csv(files[i],header=T,sep="\t",stringsAsFactors=F) 

	perm<-read.csv(paste0(dir,'perm_',names[i],'_merged') ,header=T,sep="\t")
	permP<-data.frame(SNP=perm$SNP,PermP=perm$Pvalue)
	file<-merge(file,permP,by="SNP",all.x=T) 
	file$PermP<-as.numeric(as.character(file$PermP)) 
	file$Pvalue<-as.numeric(as.character(file$Pvalue)) 
	file$TechKinPvalue<-as.numeric(as.character(file$TechKinPvalue) ) 
	gis<-getGIs(file,names[i]) 
	print(gis)
	if(i==1) inflation<-gis
	if(i>1)inflation<-data.frame(rbind(gis,inflation) ) 
}

write.table(inflation,paste0(bDir,"CaseControlResults/inflation") ,col.names=T,row.names=F,quote=F,sep="\t") 

data <- melt(inflation, id=c("Group","MAF"))
byModel<-aggregate(value~MAF+variable,data=tra,mean ) 
colnames(byModel)[grep("value",colnames(byModel))]<-"Inflation"
p <- ggplot(byModel, aes(x=MAF, y=Inflation,colour=variable)) 
p<-p + geom_line()
p<-p+ ggtitle("Genomic Inflation across UCLex Cohorts by Model")
pdf(paste0(bDir,"CaseControlResults/Model_genomic_control_lambdas.pdf") ) 
print(p)
dev.off() 




#p <- ggplot(data, aes(x=MAF, y=value,colour=variable)) 
#p + geom_line(aes(linetype = variable))
