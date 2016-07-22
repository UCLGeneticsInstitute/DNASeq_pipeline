files<-list.files('/scratch2/vyp-scratch2/cian/UCLex_February2015/KinshipDecomposition_all/',full.names=T,pattern='indi.res') 
files<-files[grep('TK_maf_0.00001_callRate_0.00001.indi.res',files)]

cohort.list<-c('Levine','Davina','Hardcastle','IoO_','IoN','IoOFFS','IoONov2013','IoOPanos','Kelsell','LambiaseSD',
'Lambiase_','LayalKC','Manchester','Nejentsev','PrionUnit','Prionb2','Shamima','Sisodiya','Syrris','Vulliamy','WebsterURMD','gosgene_')

fam<-read.table("/scratch2/vyp-scratch2/cian/UCLex_February2015/allChr_snpStats.fam",header=F,sep="\t")

for(i in 1:length(cohort.list))
{
	file<-files[grep(cohort.list[i],files)]
	tra<-read.table(file,header=T,sep="\t")

	if(i==1)
	{
		dat<-data.frame(matrix(nrow=nrow(tra),ncol=(length(cohort.list)+2)))
		dat[,1:2]<-tra[,1:2]
	}
	dat[,i+2]<-tra$Residual

}
write.table(dat,'TK_residual_pheno_file',col.names=F,row.names=F,quote=F,sep="\t")
write.table(cohort.list,'Res-pheno-file-order',col.names=F,row.names=F,quote=F,sep="\t")

