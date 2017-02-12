
################################################################################################################################
### STEP3: From "individual-separated SV-TXT files" to "merged SV-TXT files" for each SVcaller (Parallel calculating by chrs)###
################################################################################################################################

rm(list=ls());options(stringsAsFactors=FALSE)

args=commandArgs(TRUE)
i=as.numeric(args[2])
chr=args[1]

setwd("/SAN/vyplab/NCMD/b37/SV-MERGED")
SVcallers_list = read.table(file="0_SVcallers_list.txt",header=F)
Inds_list = read.table(file="0_Inds_list.txt",header=F)
max_offset = 1000 # max. breakpoint offset between two SVs #
#min_overlap = 0.8 # min. reciprocal overlap between two SVs #

SVcaller = SVcallers_list[i,1]
myfolder=paste(SVcaller,'/4_MERGED_TXT/',SVcaller,"_AllInds_Chrs",sep="")
if(paste(SVcaller,"_AllInds_Chrs",sep="") %in% dir(paste(SVcaller,'/4_MERGED_TXT/',sep=""))==FALSE) dir.create(myfolder) 

for(j in 1:dim(Inds_list)[1]){
	
	Ind = Inds_list[j,1]		
	Ind_filename = list.files(path=paste(SVcaller,"/4_MERGED_TXT",sep=""),pattern=Ind)[1]
	
	if(!is.na(Ind_filename)){
		
		TEMP_vcf_ind = read.table(file=paste(SVcaller,"/4_MERGED_TXT/",SVcaller,"_",Ind,".txt",sep=""),header=T)
		TEMP_vcf_ind = subset(TEMP_vcf_ind,TEMP_vcf_ind[,1]==chr)
		
		if(j == 1){
			# (1) Initialization from the first individual #
			TEMP_vcf_soft = TEMP_vcf_ind[,c("CHROM","POS","END","SVTYPE","GT")]
			colnames(TEMP_vcf_soft)[length(colnames(TEMP_vcf_soft))] = Ind
			TEMP_vcf_soft = cbind(paste(TEMP_vcf_soft[,"CHROM"],TEMP_vcf_soft[,"POS"],TEMP_vcf_soft[,"END"],sep="_"),TEMP_vcf_soft)
			TEMP_vcf_soft = subset(TEMP_vcf_soft,!duplicated(TEMP_vcf_soft[,1]))[,-1]
			rownames(TEMP_vcf_soft) = paste(TEMP_vcf_soft[,"CHROM"],TEMP_vcf_soft[,"POS"],TEMP_vcf_soft[,"END"],sep="_")
			cat(paste(SVcaller," : ",Ind,"'s SV-TXT file has been merged into the SV result","\n",sep=""))
		}else{
			if(i %in% c(1,6)){
				# (2) When SVcaller is Delly/Wham, directly merge the new SV genotypes with old SV file #
				TEMP_vcf_soft = cbind(TEMP_vcf_soft,TEMP_vcf_ind[,"GT"])
				colnames(TEMP_vcf_soft)[length(colnames(TEMP_vcf_soft))] = Ind
			}else{
				# (3) When SVcaller is other callers, implement different merge #
				TEMP_vcf_soft = cbind(TEMP_vcf_soft,"NA")
				colnames(TEMP_vcf_soft)[length(colnames(TEMP_vcf_soft))] = Ind
				TEMP_vcf_ind = cbind(paste(TEMP_vcf_ind[,"CHROM"],TEMP_vcf_ind[,"POS"],TEMP_vcf_ind[,"END"],sep="_"),TEMP_vcf_ind)
				TEMP_vcf_ind = subset(TEMP_vcf_ind,!duplicated(TEMP_vcf_ind[,1]))[,-1]
				rownames(TEMP_vcf_ind) = paste(TEMP_vcf_ind[,"CHROM"],TEMP_vcf_ind[,"POS"],TEMP_vcf_ind[,"END"],sep="_")
				
				for(k in 1:dim(TEMP_vcf_ind)[1]){
					
					# FILTER: max. breakpoint offset, -/+ 1000bp #
					Old_SVs = subset(TEMP_vcf_soft,TEMP_vcf_soft[,"CHROM"]==TEMP_vcf_ind[k,"CHROM"] & TEMP_vcf_soft[,"POS"] %in% ((TEMP_vcf_ind[k,"POS"]-max_offset):(TEMP_vcf_ind[k,"POS"]+max_offset)) & TEMP_vcf_soft[,"END"] %in% ((TEMP_vcf_ind[k,"END"]-max_offset):(TEMP_vcf_ind[k,"END"]+max_offset)) & TEMP_vcf_soft[,"SVTYPE"]==TEMP_vcf_ind[k,"SVTYPE"])
					
					# (3-1) There is old SVs could be merged into for the new SV (Just recode the "GT", "POS" and "END" in the NA columns) #
					if(dim(Old_SVs)[1] > 0){												
						if(dim(Old_SVs)[1] > 1){
							TEMP_vcf_soft = TEMP_vcf_soft[!rownames(TEMP_vcf_soft) %in% c(rownames(Old_SVs)[-1]),] # dim(Old_SVs)[1]>1 means the first SV-TXT file contains similar SV rows and this command just use to remove all other rows expect for the first SV row #
							Old_SVs = Old_SVs[1,] # remove the rownames of all other rows expect for the first SV row #
						}						
						TEMP_vcf_soft[rownames(Old_SVs),"POS"] = min(TEMP_vcf_soft[rownames(Old_SVs),"POS"],TEMP_vcf_ind[k,"POS"]) # Select the smaller START to use #
						TEMP_vcf_soft[rownames(Old_SVs),"END"] = max(TEMP_vcf_soft[rownames(Old_SVs),"END"],TEMP_vcf_ind[k,"END"]) # Select the bigger END to use #
						TEMP_vcf_soft[rownames(Old_SVs),Ind] = TEMP_vcf_ind[k,"GT"] # Use the GT from new SV-TXT file to fill the NA column of old SV-TXT file #
						Old_SVs_rename = paste(TEMP_vcf_soft[rownames(Old_SVs),"CHROM"],TEMP_vcf_soft[rownames(Old_SVs),"POS"],TEMP_vcf_soft[rownames(Old_SVs),"END"],sep="_")
						if(Old_SVs_rename %in% rownames(TEMP_vcf_soft)){
							TEMP_vcf_soft = TEMP_vcf_soft[!rownames(TEMP_vcf_soft) %in% c(Old_SVs_rename,rownames(Old_SVs)),] # In case for the different SVTYPE but with similar START and END #
						}else{
							rownames(TEMP_vcf_soft)[rownames(TEMP_vcf_soft)==rownames(Old_SVs)] = Old_SVs_rename # Rename that merged row of old SV-TXT file #
						}
					}else{
						# (3-2) There isn't old SVs could be merged into for the new SV (Add a new SV row in the old SV file)  #
						New_SVs = cbind(TEMP_vcf_ind[k,c("CHROM","POS","END","SVTYPE")],t(rep("NA",c(dim(TEMP_vcf_soft)[2]-5))),TEMP_vcf_ind[k,"GT"])
						colnames(New_SVs) = colnames(TEMP_vcf_soft)
						TEMP_vcf_soft = rbind(TEMP_vcf_soft,New_SVs)
					}
				}
			}
			cat(paste(SVcaller," : ",Ind,"'s SV-TXT file has been merged into the SV result","\n",sep=""))
		}
		rm(list=ls(pattern="tmp_vcf*"))
	}else{
		TEMP_vcf_soft = cbind(TEMP_vcf_soft,"NA")
		colnames(TEMP_vcf_soft)[length(colnames(TEMP_vcf_soft))] = Ind
		cat(paste(SVcaller," : ",Ind,"'s SV-TXT file is NA","\n",sep=""))
	}
	
	TEMP_vcf_soft = TEMP_vcf_soft[order(TEMP_vcf_soft[,1],TEMP_vcf_soft[,2]),]
}

write.table(TEMP_vcf_soft,file=paste(SVcaller,"/4_MERGED_TXT/",SVcaller,"_AllInds_Chrs/",SVcaller,"_AllInds_Chr",chr,".txt",sep=""),row.names=F,col.names=T,quote=F)
cat(paste(SVcaller,"############################################################","\n",sep=""))


### Merge the chr-separated SV-TXT files together for each SVcaller ###

if(length(dir(paste(SVcaller,"/4_MERGED_TXT/",SVcaller,"_AllInds_Chrs/",sep="")))==24){
	tmp_chr_all = c()
	for(chr in c(1:22,"X","Y")){
		tmp_chr = read.table(file=paste(SVcaller,"/4_MERGED_TXT/",SVcaller,"_AllInds_Chrs/",SVcaller,"_AllInds_Chr",chr,".txt",sep=""),header=T)
		tmp_chr_all = rbind(tmp_chr_all,tmp_chr)
		print(chr)
	}
	tmp_chr_all = tmp_chr_all[order(tmp_chr_all[,1],tmp_chr_all[,2]),]
	
	write.table(tmp_chr_all,file=paste(SVcaller,"/4_MERGED_TXT/",SVcaller,"_AllInds.txt",sep=""),row.names=F,col.names=T,quote=F)
	cat(paste(SVcaller,"############################################################","\n",sep=""))
}else{
	print(paste("There are ",abs(length(dir(paste(SVcaller,"/4_MERGED_TXT/",SVcaller,"_AllInds_Chrs/",sep="")))-24),"chr-separated files hasn't been generated",sep=""))
}
