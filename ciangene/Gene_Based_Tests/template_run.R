rootODir<-'/SAN/vyplab/UCLex/mainset_July2016/cian/'
release<-'July2016'

pheno.matching<-data.frame(read.csv('/SAN/vyplab/UCLex/scripts/DNASeq_pipeline/ciangene/Support/phenotype.info.csv',header=FALSE))

oDir<-paste0(rootODir,'SKATnew/')
if(!file.exists(oDir)) dir.create(oDir)

#for(phen in 73:nrow(pheno.matching))
for(phen in 73:110)
{
	script.out<-paste0(oDir,pheno.matching[phen,1],'.R')
	oData<-paste0('load("',paste0(oDir,pheno.matching[phen,1],'.RData"'),')')
	write.table(oData,script.out,col.names=F,row.names=F,quote=F,sep='\t')	
	oDirr<-paste0('oDir<-"',oDir,'"')
	write.table(oDirr,script.out,col.names=F,row.names=F,quote=F,sep='\t',append=TRUE)	
	runR='sh /cluster/project8/vyp/cian/scripts/bash/runRonCluster.sh'
	file.append( script.out, 'SKAT_uclex.template.R') 
	system(paste(runR,script.out))
}
