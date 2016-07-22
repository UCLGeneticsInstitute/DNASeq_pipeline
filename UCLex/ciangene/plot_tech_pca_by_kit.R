library(ggplot2) 
library(gridExtra)
library(grid)
release<-"July2015"
bDir<-paste0('/scratch2/vyp-scratch2/cian/UCLex_',release,"/") 
#source("multiplot.R") 
file<-read.table("/cluster/project8/vyp/cian/data/UCLex/UCLexInfo/uclex-samples.csv",header=T,sep=",")

prepData<-function(pcs)
{
	pcs$Kit<-NA
	kits<-unique(file$sequencing.platform) 
	for( i in 1:length(kits) ) 
	{
		samples<-file$sample[which(file$sequencing.platform == kits[i])]
		pcs$Kit[pcs[,1]%in%samples]<-kits[i]
	}

	seq<-read.table("/cluster/project8/vyp/cian/data/UCLex/UCLexInfo/Some_Sample_Capture_Kit_Info",header=T,sep="\t")
	seq$enrichment.kit[grep("unknown",seq$enrichment.kit)]<-NA

	kits<-unique(seq$enrichment.kit)
	kits<-kits[!is.na(kits) ]
	for( i in 1:length(kits) ) 
	{
		samples<-seq[,1][which(seq$enrichment.kit == kits[i])]
		pcs$Kit[pcs[,1]%in%samples ]<-kits[i]
	}

	data<-data.frame(PC1=pcs$V3,PC2=pcs$V4,Kit=pcs$Kit) 
	return(data) 
}
	 
techPCs<-read.table(paste0(bDir,"TechPCs.vect"),header=F)
depthPCs<-read.table(paste0(bDir,"/DepthPCs.vect"),header=F)
tech.short<-data.frame(techPCs[,1:4])
depth.short<-data.frame(depthPCs[,1:4])
colnames(tech.short)<-c('Sample1','Sample2','PC1','PC2') 
colnames(depth.short)<-c('Sample1','Sample2','PC1','PC2') 

techPCA<-merge(tech.short,file,by.x='Sample1',by.y='sample')
depthPCA<-merge(depth.short,file,by.x='Sample1',by.y='sample')


grid_arrange_shared_legend <- function(...) {
    plots <- list(...)
    g <- ggplotGrob(plots[[1]] + theme(legend.position="bottom"))$grobs
    legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
    lheight <- sum(legend$height)
    grid.arrange(
        do.call(arrangeGrob, lapply(plots, function(x)
            x + theme(legend.position="none") + theme(legend.text=element_text(size=.7)) )),
        legend,
        ncol = 1,
        heights = unit.c(unit(1, "npc") - lheight, lheight))
}

pdf("test.pdf")
grid_arrange_shared_legend(top.plot ,bottom.plot)
dev.off()


exit
dat<-data.frame(sample=techPCA[,1],kit=techPCA$Kit) 
tra<-merge(file,dat,by="sample") 



#techPCA<-prepData(techPCs) 
#depthPCA<-prepData(depthPCs) 

old<-FALSE
if(old)
{
tech<-qplot(PC1,PC2,colour=sequencing.platform,data=techPCA,main="Missingness PCA")
depth<-qplot(PC1,PC2,colour=sequencing.platform,data=depthPCA,main="Depth PCA")

pdf(paste0(bDir,"identifying_capture_tech_by_pca.pdf"), onefile = TRUE)
#print(multiplot(tech,depth,col=1) ) 
top.plot <- qplot(PC1,PC2,colour=sequencing.platform,data=techPCA,main="Missingness PCA") 
bottom.plot <- qplot(PC1,PC2,colour=sequencing.platform,data=depthPCA,main="Depth PCA")
grid.arrange(top.plot, bottom.plot,ncol=3) #layout_matrix=rbind(c(1,2))) 
dev.off()
}

