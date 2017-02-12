
############################################################################################################################
### STEP4: From each SVcaller's merged SV-TXT files for the whole merged SV-TXT (Parallel calculating by chrs)###
############################################################################################################################

rm(list=ls());options(stringsAsFactors=FALSE)

args=commandArgs(TRUE)
chr=args[1]

setwd("/SAN/vyplab/NCMD/b37/SV-MERGED")
SVcallers_list = read.table(file="0_SVcallers_list.txt",header=F)
Inds_list = read.table(file="0_Inds_list.txt",header=F)
max_offset = 1000 # max. breakpoint offset between two SVs #
#min_overlap = 0.8 # min. reciprocal overlap between two SVs #

myfolder="Allcallers_AllInds_Chrs"
if(myfolder %in% dir()==FALSE) dir.create(myfolder) 

for(i in 1:dim(SVcallers_list)[1]){

	### STEP1. Load the merged VCF file of each SVcaller ###
	TEMP_vcf_soft = read.table(file=paste(SVcallers_list[i,1],"/4_MERGED_TXT/",SVcallers_list[i,1],"_AllInds.txt",sep=""),header=T)		
	TEMP_vcf_soft = cbind(paste(TEMP_vcf_soft[,"CHROM"],TEMP_vcf_soft[,"POS"],TEMP_vcf_soft[,"END"],sep="_"),TEMP_vcf_soft)
	TEMP_vcf_soft = subset(TEMP_vcf_soft,!duplicated(TEMP_vcf_soft[,1]))[,-1]
	rownames(TEMP_vcf_soft) = paste(TEMP_vcf_soft[,"CHROM"],TEMP_vcf_soft[,"POS"],TEMP_vcf_soft[,"END"],sep="_")
	colnames(TEMP_vcf_soft)[5:dim(TEMP_vcf_soft)[2]] = paste(SVcallers_list[i,1],colnames(TEMP_vcf_soft)[5:dim(TEMP_vcf_soft)[2]],sep="_")	
	TEMP_vcf_soft = subset(TEMP_vcf_soft,TEMP_vcf_soft[,"CHROM"]==chr)
	
	### STEP2. Start to merge all the VCF files together ###
	if(i == 1){
		# (1) Initialization from the first SVcaller #
		TEMP_vcf_soft_all = TEMP_vcf_soft
	}else{
		# (2) When SVcaller is other callers, implement different merge #
		TEMP_vcf_soft_all = cbind(TEMP_vcf_soft_all,t(rep("NA",dim(Inds_list)[1])))
		Add_cols = c(((dim(Inds_list)[1]+4)+(dim(Inds_list)[1]*(i-2)+1)):((dim(Inds_list)[1]+4)+dim(Inds_list)[1]*(i-1)))
		Source_cols = c(5:dim(TEMP_vcf_soft)[2])
		colnames(TEMP_vcf_soft_all)[Add_cols] = colnames(TEMP_vcf_soft)[Source_cols]
		
		for(j in 1:dim(TEMP_vcf_soft)[1]){
			
			# FILTER: max. breakpoint offset, -/+ 1000bp #
			Old_SVs = subset(TEMP_vcf_soft_all,TEMP_vcf_soft_all[,"CHROM"]==TEMP_vcf_soft[j,"CHROM"] & TEMP_vcf_soft_all[,"POS"] %in% ((TEMP_vcf_soft[j,"POS"]-max_offset):(TEMP_vcf_soft[j,"POS"]+max_offset)) & TEMP_vcf_soft_all[,"END"] %in% ((TEMP_vcf_soft[j,"END"]-max_offset):(TEMP_vcf_soft[j,"END"]+max_offset)) & TEMP_vcf_soft_all[,"SVTYPE"]==TEMP_vcf_soft[j,"SVTYPE"])
			
			# (2-1) There is old SVs could be merged into for the new SV (Just recode the "GT", "POS" and "END") #
			if(dim(Old_SVs)[1]>0){
				if(dim(Old_SVs)[1] > 1){
					TEMP_vcf_soft_all = TEMP_vcf_soft_all[!rownames(TEMP_vcf_soft_all) %in% c(rownames(Old_SVs)[-1]),] # dim(Old_SVs)[1]>1 means the first SV-TXT file contains similar SV rows and this command just use to remove all other rows expect for the first SV row #
					Old_SVs = Old_SVs[1,] # remove the rownames of all other rows expect for the first SV row #
				}	
				TEMP_vcf_soft_all[rownames(Old_SVs),"POS"] = round(mean(c(TEMP_vcf_soft_all[rownames(Old_SVs),"POS"],TEMP_vcf_soft[j,"POS"])))
				TEMP_vcf_soft_all[rownames(Old_SVs),"END"] = round(mean(c(TEMP_vcf_soft_all[rownames(Old_SVs),"END"],TEMP_vcf_soft[j,"END"])))
				TEMP_vcf_soft_all[rownames(Old_SVs),Add_cols] = TEMP_vcf_soft[j,Source_cols]
			
				Old_SVs_rename = paste(TEMP_vcf_soft_all[rownames(Old_SVs),"CHROM"],TEMP_vcf_soft_all[rownames(Old_SVs),"POS"],TEMP_vcf_soft_all[rownames(Old_SVs),"END"],sep="_")
				if(Old_SVs_rename %in% rownames(TEMP_vcf_soft_all)){
					TEMP_vcf_soft_all = TEMP_vcf_soft_all[!rownames(TEMP_vcf_soft_all) %in% c(Old_SVs_rename,rownames(Old_SVs)),] # In case for the different SVTYPE but with similar START and END #
				}else{
					rownames(TEMP_vcf_soft_all)[rownames(TEMP_vcf_soft_all)==rownames(Old_SVs)] = Old_SVs_rename # Rename that merged row of old SV-TXT file #
				}
			
			}else{
				# (2-2) There isn't old SVs could be merged into for the new SV (Add a new SV row in the old SV file)  #
				New_SVs = cbind(TEMP_vcf_soft[j,c("CHROM","POS","END","SVTYPE")],t(rep("NA",(i-1)*dim(Inds_list)[1])),TEMP_vcf_soft[j,Source_cols])
				colnames(New_SVs) = colnames(TEMP_vcf_soft_all)
				TEMP_vcf_soft_all = rbind(TEMP_vcf_soft_all,New_SVs)
			}
			cat(paste(SVcallers_list[i,1],"-",j," ",sep=""),sep="")
		}
	}
	
	cat(paste("\n",SVcallers_list[i,1],"'s VCF file has been merged into the final SV result","\n",sep=""))	
}

TEMP_vcf_soft_all = TEMP_vcf_soft_all[order(TEMP_vcf_soft_all[,1],TEMP_vcf_soft_all[,2]),]
NEW_Ind_order = c("CHROM","POS","END","SVTYPE",paste(rep(SVcallers_list[,1],dim(Inds_list)[1]),rep(Inds_list[,1],each=dim(SVcallers_list)[1]),sep="_"))
TEMP_vcf_soft_all = TEMP_vcf_soft_all[,NEW_Ind_order]
write.table(TEMP_vcf_soft_all,file=paste("Allcallers_AllInds_Chrs/Allcallers_AllInds_Chr",chr,".txt",sep=""),row.names=F,col.names=T,quote=F)
cat(paste("All SVcaller at chr ",chr," have been merged together! ^_^","\n",sep=""))


### Merge the chr-separated SV-TXT files together for each SVcaller ###

if(length(dir("Allcallers_AllInds_Chrs/"))==24){
	tmp_chr_all = c()
	for(chr in c(1:22,"X","Y")){
		tmp_chr = read.table(file=paste("Allcallers_AllInds_Chrs/Allcallers_AllInds_Chr",chr,".txt",sep=""),header=T)
		tmp_chr_all = rbind(tmp_chr_all,tmp_chr)
		print(chr)
	}
	tmp_chr_all = tmp_chr_all[order(tmp_chr_all[,1],tmp_chr_all[,2]),]
	
	write.table(tmp_chr_all,file="Allcallers_AllInds.txt",row.names=F,col.names=T,quote=F)
	cat(paste("All SVcallers has been merged ############################################################","\n",sep=""))
}else{
	print(paste("There are ",abs(length(dir("Allcallers_AllInds_Chrs/"))-24),"chr-separated files hasn't been generated",sep=""))
}


### Calculate the agreement table ###

Agreement_table = as.data.frame(matrix(NA,dim(tmp_chr_all)[1],4+dim(Inds_list)[1]))
Agreement_table[,1:4] = tmp_chr_all[,1:4]
colnames(Agreement_table) = c("CHROM","POS","END","SVTYPE",t(Inds_list[,1]))

for(k in 5:(dim(Agreement_table)[2])){
	tmp = tmp_chr_all[,(4+dim(SVcallers_list)[1]*(k-4)-c(dim(SVcallers_list)[1]-1)):(4+dim(SVcallers_list)[1]*(k-4))]	
	Agreement_table[,k] = apply(tmp,1,function(x){return(length(x[!is.na(x)])/dim(SVcallers_list)[1])})	
	print(k)
}
write.table(Agreement_table,file=paste("Allcallers_AllInds_Agreement.txt",sep=""),row.names=F,col.names=T,quote=F)
