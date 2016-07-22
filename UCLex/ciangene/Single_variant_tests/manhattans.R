source("/cluster/project8/vyp/cian/scripts/r/manhattan.R") 
bDir<-paste0("/scratch2/vyp-scratch2/cian/UCLex_",release,'/')
dir<-paste0(bDir,"FastLMM_Single_Variant_all_phenos/") 
files<-list.files(dir,pattern="final",full.names=T) 
names<-gsub(basename(files),pattern="_.*",replacement="") 

#pdf("Manhattans.pdf") 
for( i in 1:length(names))
{
	print(names[i]) 
	file<-read.csv(files[i],header=T,sep="\t",stringsAsFactors=F) 
	tra<-subset(file,file$ExtCtrl_MAF>0.01) 
	tra$TechKinPvalue<-as.numeric(as.character(tra$TechKinPvalue)) 
	dat<-data.frame(CHR=tra$Chr,BP=tra$Start,P=tra$TechKinPvalue) 
	dat$CHR[grep("X",dat$CHR)]<-23
	dat$CHR<-as.numeric(dat$CHR) 
	
	oFile<-paste0(bDir,names[i],"_manhattan.png") 	
	png(oFile) 
	manhattan(dat,main=names[i]) 
	dev.off()
}
#dev.off()


