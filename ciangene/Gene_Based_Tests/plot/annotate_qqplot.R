source("qqplot.R") ## modified snpStats qq.chisq function to return obsvd and expctd values. 

labelQQplot <- function(pval.labels, qqplotOut ,pvals.raw, nb.pvals=5, col='red', GenesToLabel=3) 
{
	dat <- data.frame(pval.labels, pvals.raw) 
	qqplotOut<-qqplotOut[-c(1:3)]
	obsvd <- qqplotOut[1: (length(qqplotOut)/2)]
	expctd<- qqplotOut[((length(qqplotOut)/2)+1): length(qqplotOut)]
	plotty <- data.frame(obsvd,expctd) 
	tra <- plotty[order(-plotty$expctd),]
	plot.data <- data.frame(labels=pval.labels,pvals=pvals.raw) 
	plot.data<-plot.data[order(plot.data$pvals),]

	if(is.numeric(GenesToLabel)) 
	{
		message("Labelling top genes")
		plot.vals <- data.frame(plot.data,tra)[1:nb.pvals,]
	} else
	{
		message("Labelling specified genes") 
		plot.vals.tra<-data.frame(plot.data,tra)
		plot.vals<-plot.vals.tra[plot.vals.tra$labels%in%genes,]
	}

	upper.labs <- plot.vals[seq(from=1,by=2, to=nrow(plot.vals)),]
	lower.labs <- plot.vals[seq(from=2,by=2, to=nrow(plot.vals)),]
	
	x.offset=1
	y.offset=0
	text(upper.labs$expctd-x.offset,upper.labs$obsvd+y.offset, labels=upper.labs$labels, col = col, adj=c(1,0.5),cex=0.6 ) 
	segments(upper.labs$expctd, upper.labs$obsvd, upper.labs$expctd-x.offset, upper.labs$obsvd+y.offset, col =col) 

	text(lower.labs$expctd+x.offset,lower.labs$obsvd-y.offset, labels=lower.labs$labels, col = col, adj=c(0,0.5) ,cex=0.6) 
	segments(lower.labs$expctd, lower.labs$obsvd, lower.labs$expctd+x.offset, lower.labs$obsvd-y.offset, col =col) 
}



#file <- read.table("regressALL", header=T) 
#tst <- qq.chisq(-2*log(file$LRT_P_Perm), df=2, x.max=30)
#labelQQplot(file$Gene_name, tst, file$LRT_P_Perm) 
# genes <- c("NUP210L", "PRPF4B", "KIF9") 
