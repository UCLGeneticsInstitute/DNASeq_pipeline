getArgs <- function() {
  myargs.list <- strsplit(grep("=",gsub("--","",commandArgs()),value=TRUE),"=")
  myargs <- lapply(myargs.list,function(x) x[2] )
  names(myargs) <- lapply(myargs.list,function(x) x[1])
  return (myargs)
}

#release <- 'June2016'
#rootODir<-'/cluster/project8/vyp/cian/data/UCLex'

myArgs <- getArgs()

if ('rootODir' %in% names(myArgs))  rootODir <- myArgs[[ "rootODir" ]]
if ('release' %in% names(myArgs))  release <- myArgs[[ "release" ]]

######################
library(plotly)


oDir <- paste0(rootODir, "/")
print(oDir)
file<-read.table(paste0(oDir,"UCLex",release,"_OneKG_merged_pca.vect"),header=F,sep=" ")
#file<-read.table(paste0(oDir,"UCLex_OneKG_merged_outliers_removed_pca.vect"),header=F,sep=" ")

sample.inf<-read.table("/SAN/vyplab/UCLex/data/OneKG_sample_ids.tab",header=T,sep="\t")
populations<-unique(sample.inf$pop)
pcs<-data.frame(file[,1:7])
pcs$Population<-"UCLex"
pcs$Caucasian<-FALSE
caucasian.groups<-c("CEU","TSI","GBR","IBS","FIN")
onekg.data<-data.frame(matrix(nrow=length(populations),ncol=3))
colnames(onekg.data)<-c("Population","PC1centre","PC2centre")
onekg.data[,1]<-populations

plot.oneKG<-TRUE
plot.uclex<-TRUE

png(paste0(oDir,"plots/Ancestry_pca.png"))

for(i in 1:length(populations))
{
	pop<-populations[i]
	pop.samples<-sample.inf$namess[grep(pop,sample.inf$pop)]
	pop.data<-pcs[pcs[,1]%in%pop.samples,]
	message(paste("There are", length(pop.samples),"in",pop))
	pcs$Population[pcs[,1]%in%pop.samples]<-as.character(pop)
	mean.pc.one<-mean(pop.data[,3])
	mean.pc.two<-mean(pop.data[,4])
	onekg.data[i,2]<-mean.pc.one
	onekg.data[i,3]<-mean.pc.two

	if(i==1)
	{
		x.min<-min(pcs[,3])
		x.max<-max(pcs[,3])
		y.min<-min(pcs[,4])
		y.max<-max(pcs[,4])
		if(plot.oneKG)col<-"yellow" else col<-"white"

		plot(pop.data[,3],pop.data[,4],
				xlim=c(x.min,x.max),
				ylim=c(y.min,y.max),
				xlab="PC1",
				ylab="PC2",
				main=paste(paste0("UCLex",release), "1000G PCA"),
				#main = ('UCL-ex 1000g PCA'),
				col=col
				)

	} else
	{
		points(pop.data[,3],pop.data[,4],col=col)
	}
	#text(mean.pc.one,mean.pc.two,pop,col="red")
}
uclex.samples<-pcs[pcs$Population=="UCLex",]
ucl.colour<-"dodgerblue4"
if(plot.uclex)
{
points(uclex.samples[,3],uclex.samples[,4],pch=4,col=ucl.colour)
}
#text(onekg.data$PC1centre,onekg.data$PC2centre,pop.data$Population,col="red")

## Redo this bit so i get text labels on top
labels<-data.frame(matrix(nrow=length(populations),ncol=3))
labels[,1]<-populations
for(i in 1:length(populations))
{
	pop<-populations[i]
	pop.samples<-sample.inf$namess[grep(pop,sample.inf$pop)]
	pop.data<-pcs[pcs[,1]%in%pop.samples,]
	pcs$Population[pcs[,1]%in%pop.samples]<-as.character(pop)
	mean.pc.one<-mean(pop.data[,3])
	mean.pc.two<-mean(pop.data[,4])
	labels[i,2]<-mean.pc.one
	labels[i,3]<-mean.pc.two
#text(mean.pc.one,mean.pc.two,pop,col="red")
}
write.table(labels,paste0(oDir,'pca.ancestry.labels.csv'),col.names=F,row.names=F,quote=F,sep='\t')



pcs$Caucasian[pcs$Population%in%caucasian.groups]<-TRUE
cauc.pc1.min<-min(pcs[pcs$Caucasian,3])*1.5
cauc.pc1.max<-max(pcs[pcs$Caucasian,3])*.3
cauc.pc2.min<-min(pcs[pcs$Caucasian,4])
cauc.pc2.max<-max(pcs[pcs$Caucasian,4])*1.4
segments(cauc.pc1.min,cauc.pc2.min,cauc.pc1.min,cauc.pc2.max,col="forestgreen") # bottom left to top left
segments(cauc.pc1.min,cauc.pc2.min,cauc.pc1.max,cauc.pc2.min,col="forestgreen")# bottom left to bottom right
segments(cauc.pc1.max,cauc.pc2.min,cauc.pc1.max,cauc.pc2.max,col="forestgreen")# bottom right to top right
segments(cauc.pc1.min,cauc.pc2.max,cauc.pc1.max,cauc.pc2.max,col="forestgreen")#

legend('topleft',legend=c("OneKG","UCLex"),col=c(col,ucl.colour),pch=16,pt.bg=c(col,ucl.colour)) 
dev.off()

for(i in 1:nrow(uclex.samples))
{
	x<-uclex.samples$V3[i]
	y<-uclex.samples$V4[i]
	if(x>cauc.pc1.min&x<cauc.pc1.max&y>cauc.pc2.min&y<cauc.pc2.max)uclex.samples$Caucasian[i]<-TRUE
}

write.table(uclex.samples,paste0(oDir,"UCLex_samples_ancestry"),col.names=T,row.names=F,quote=F,sep="\t")




pca.plot<-ggplot(pcs, aes(V3, V4, colour=Population) ) +geom_point(shape=1) +  ggtitle(paste(date(),'UCLex ancestry PCA'))
pca.plot <- pca.plot +geom_rect(data=pcs, mapping=aes(xmin=cauc.pc1.min, xmax=cauc.pc1.max, ymin=cauc.pc2.min, ymax=cauc.pc2.max),alpha=0,colour='black') 
pca.plot+annotate("text", x = labels[,2], y = labels[,3], label = labels[,1]) 
g<-ggplotly(pca.plot)
p <- plotly_build(g) 
p$filename <- 'UCLex_ancestry_pca'
r <- plotly_POST(p)
knit_print.plotly(r, options=list())


##  batches<-t(data.frame(strsplit(pcs[,1],'_')) ) - maybe finish this to have different research groups represented on graph as separate shapes. 


pcs$Group<-NA
colnames(pcs)[3:4]<-c('PC1','PC2')
pca.plot<-plot_ly(pcs, x = PC1, y = PC2, text = V1, mode="markers",color=Population, main='UCLex 1000G PCA')#+geom_point(shape=1) #+  ggtitle(paste(date(),'UCLex ancestry PCA'))
g<-ggplotly(pca.plot)
p <- plotly_build(pca.plot) 
p$filename <- 'UCLex_ancestry_pca_names'
r <- plotly_POST(p)
knit_print.plotly(r, options=list())

htmlwidgets::saveWidget(as.widget(p), paste0(oDir,"graph.html")) 
htmlwidgets::saveWidget(as.widget(p), "graph2.html")

save.image(file=paste0(oDir,'pca.plot.RData')) 
save.image(file='pca.plot.RData')



#### The important stuff. 
pca.plot<-plot_ly(pcs, x = PC1, y = PC2, text = V1, mode="markers",color=Population, title='UCLex 1000G PCA')%>%
  layout(title = paste(date(),'UCLex ancestry PCA')) 

p <- plotly_build(pca.plot) 

p <- layout(p,
             shapes = list(
               list(type = "rect",
                 fillcolor = "blue", line = list(color = "blue"), opacity = 0.2,
                 x0 = cauc.pc1.min, x1 = cauc.pc1.max, xref = "x",
                 y0 = cauc.pc2.min, y1 = cauc.pc2.max, yref = "y")))

htmlwidgets::saveWidget(as.widget(p), "UCLex_ancestry_pca.html")

