suppressPackageStartupMessages(library("optparse"))

option_list <- list(
	make_option(c("--depthDir"), default=NULL,help="Chromosome"),
 	make_option(c("--bamDir"), default='/SAN/vyplab/UCLex_raw/',help="Chromosome")
 	)
opt <- parse_args(OptionParser(option_list=option_list))

message("Starting Argument checks...\n", stderr())

depthDir<-opt$depthDir
bamDir<-opt$bamDir

files<-list.files(depthDir,full.names=T,pattern='txt$')
if(length(which(grepl('gendered',files)))>0) files<-files[-grep('gendered',files) ]
print(files)

x<-data.frame(matrix(nrow=length(files),ncol=2))
x[,1]<-files

for(i in 1:length(files))
{
    dat<-read.table(files[i])
    dat$Ratio<-dat[,3]/dat[,4]
    print(dat$Ratio[dat[,1]=='Y'])
    x[i,2]<-dat$Ratio[dat[,1]=='Y']
}
x$MALE<-TRUE
x$MALE[x[,2]<10]<-FALSE
print(table(x$MALE)) 
x$name<-file.path(bamDir,gsub(gsub(basename(x[,1]),pattern='ChrDepth/',replacement=''),pattern='.chromdepths.txt',replacement='') )
write.table(x,file.path(depthDir,'BAMS.gendered.tab'),col.names=T,row.names=F,quote=F,sep='\t')

x.out<-file.path(depthDir,'X')
if(!file.exists(x.out))  dir.create(x.out)
write.table(x$name[x$MALE],file.path(x.out,'Bams_male'),col.names=F,row.names=F,quote=F,sep='\t')
write.table(x$name[!x$MALE],file.path(x.out,'Bams_female'),col.names=F,row.names=F,quote=F,sep='\t')
