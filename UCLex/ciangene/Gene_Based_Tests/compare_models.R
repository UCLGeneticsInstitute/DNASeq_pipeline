library(snpStats) 
cohort<-"IoN"
dDir<-'/scratch2/vyp-scratch2/cian/UCLex_February2015/LDAK_gene_tests_all_phenos/'

files <- list.files(dDir, pattern = cohort) 

inBase <- paste0(dDir, "no_kin_", cohort, "_.0000001_.00000001/regressALL") 
base <- read.table(inBase, header=T) 

tech <- read.table("/scratch2/vyp-scratch2/cian/UCLex_February2015/LDAK_gene_tests_all_phenos//with_kin_IoN_.0000001_.00000001/regressALL", header=T) 
tech.small <- data.frame(tech$Gene_name, tech$LRT_P_Perm) 
colnames(tech.small)<- gsub(colnames(tech.small), pattern = "tech\\.", replacement ="")
colnames(tech.small)[2]<-"TechKinPvalue"
res <- read.table("/scratch2/vyp-scratch2/cian/UCLex_February2015/LDAK_gene_tests_all_phenos/no_kin_IoN_.0000001_.00000001_res/regressALL", header=T) 
res.small <- data.frame(res$Gene_name, res$LRT_P_Perm) 
colnames(res.small)<- gsub(colnames(res.small), pattern = "res\\.", replacement ="")
colnames(res.small)[2]<-"ResKinPvalue"
base.tech <- merge(base, tech.small, by="Gene_name") 
base.tech.res <- merge(base.tech, res.small, by="Gene_name") 

png(paste0(dDir,"model_comparison.png") ) 
par(mfrow=c(2,2))  
qq.chisq(-2*log(base.tech.res$LRT_P_Perm), df=2, x.max=30, main = "IoN Base 10% missingness") 
qq.chisq(-2*log(base.tech.res$TechKinPvalue), df=2, x.max=30, main = "IoN TechKin 10% missingness") 
qq.chisq(-2*log(base.tech.res$ResKinPvalue), df=2, x.max=30, main = "IoN Res 10% missingness") 
dev.off()


