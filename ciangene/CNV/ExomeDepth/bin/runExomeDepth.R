Rbin='/usr/local/bin/R'
bamDir<-"/SAN/vyplab/UCLex/BAM/"
bams<-list.files(bamDir,pattern="bam$",full.names=T,recursive=T) 
bams<-bams[grep('IRDC',bam)]
runSh='sh /cluster/project8/vyp/cian/scripts/bash/runBashCluster.sh'
batches<-basename(dirname(dirname(bams))) 
uniq.batches<-unique(batches)

bDir<-'/cluster/project8/vyp/cian/data/Cardio/HCM/CNV/ExomeDepth/'
outDir<-paste0(bDir,"Calls/") 
if(!file.exists(outDir))dir.create(outDir) 
max.group.size<-12

finish<-function(out,samples,name)
{
	if(length(which(is.na(samples)))>0) stop("There are na samples") 
	sample.list<-paste0(out,"_samples")
	write.table(samples,sample.list,col.names=F,row.names=F,quote=F,sep="\t")
	sample.list<-paste0(basename(out),"_samples") 
	run<-paste(Rbin," CMD BATCH --no-save --no-restore",paste0("--samples=",outDir,basename(out),"_samples"), 	paste0("--Out=",name), 
			paste0(bDir,"exomeDepth_rscript.R") ,paste0(outDir,basename(out),".Rout")  
			) 
	write.table(run,out,col.names=F,row.names=F,quote=F,sep="\t") 	
	system(paste(runSh,out)) 
}
for(batch in 1:length(uniq.batches))
{
	samples<-bams[grep(uniq.batches[batch],batches)]
	if(length(samples)>max.group.size) 
	{
		nb.groups<-ceiling(length(samples)/max.group.size) 
		for(group in 1:nb.groups) 
		{
			samples<-samples[!is.na(samples)]
			sample.set<-samples[1:max.group.size] 
	
			if(length(samples)-max.group.size*2<0)
			{	
				sample.set<-samples[!is.na(samples)]
				print(paste0("Making a leftover group of ", length(samples) ) ) 
				out<-paste0("./Calls/",batches[grep(uniq.batches[batch],batches)][1],"_subgroup_", group) 
				name<-paste0(batches[grep(uniq.batches[batch],batches)][1],"_subgroup_", group) 
				finish(out,sample.set,name)
				break
			}
			print(paste("There are",length(sample.set),"samples")) 
			samples[samples%in%sample.set]<-NA
			out<-paste0("./Calls/",batches[grep(uniq.batches[batch],batches)][1], "_subgroup_",group) 
			name<-paste0(batches[grep(uniq.batches[batch],batches)][1], "_subgroup_",group) 
			finish(out,sample.set,name)
		}
	}else 
	{
		out<-paste0("./Calls/",batches[grep(uniq.batches[batch],batches)][1]) 
		name<-paste0(batches[grep(uniq.batches[batch],batches)][1]) 
		finish(out,samples,name)
	}
} 
