dat<-read.csv("/scratch2/vyp-scratch2/cian/UCLex_July2015/FastLMM_Single_Variant_all_phenos/Lambiase_single_variant_vs_UCLex.csv") 
write.table(dat$SNP,"Lambiase_snps",col.names=F,row.names=F,quote=F,sep="\t") 
source("check_Lambiase.R") 

