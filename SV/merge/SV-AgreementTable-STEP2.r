
###############################################################################################
### STEP2: From "individual-separated VCF/BCF files" to "individual-separated SV-TXT files" ###
###############################################################################################

rm(list=ls());options(stringsAsFactors=FALSE)

args=commandArgs(TRUE)
i=as.numeric(args[1])

setwd("/SAN/vyplab/NCMD/b37/SV-MERGED")
SVcallers_list = read.table(file="0_SVcallers_list.txt",header=F)
Inds_list = read.table(file="0_Inds_list.txt",header=F)

SVcaller = SVcallers_list[i,1]

for(j in 1:dim(Inds_list)[1]){
	
	Ind = Inds_list[j,1]		
	Ind_filename = list.files(path=paste(SVcaller,"/2_SEPARATED",sep=""),pattern=Ind)[1]
	
	if(!is.na(Ind_filename)){
		system(paste("bcftools view -h ",SVcaller,"/2_SEPARATED/",Ind_filename," | tail -1 | sed 's/#//g' > ",SVcaller,"/tmp_vcfheader.txt",sep=""))
		system(paste("bcftools view -H ",SVcaller,"/2_SEPARATED/",Ind_filename," > ",SVcaller,"/tmp_vcfbody.txt",sep=""))
		tmp_vcf = read.table(file=paste(SVcaller,"/tmp_vcfbody.txt",sep=""),header=F)
		colnames(tmp_vcf) = read.table(file=paste(SVcaller,"/tmp_vcfheader.txt",sep=""),header=F)
		
		if(SVcaller=="Delly"){
			tmp_vcf_filter = tmp_vcf
			tmp_vcf_1 = tmp_vcf_filter[,c("CHROM","POS","REF","ALT")]
			tmp_vcf_2 = do.call("rbind",strsplit(tmp_vcf_filter[,8],";"))[,c(5,2,10,11)]
			tmp_vcf_3 = do.call("rbind",strsplit(tmp_vcf_filter[,10],":"))[,1:12]
			tmp_vcf_4 = cbind(tmp_vcf_1,tmp_vcf_2,tmp_vcf_3)[,c(1,2,5,3,4,6:20)]
			TEMP_vcf_ind = cbind(tmp_vcf_4[,1:6],"NA",tmp_vcf_4[,7:20])
			colnames(TEMP_vcf_ind) = c("CHROM","POS","END","REF","ALT","SVTYPE","SVLEN","CIPOS","CIEND","GT","GL","GQ","FT","RCL","RC","RCR","CN","DR","DV","RR","RV")
			TEMP_vcf_ind[,3] = gsub("END=","",TEMP_vcf_ind[,3])
			TEMP_vcf_ind[,6] = gsub("SVTYPE=","",TEMP_vcf_ind[,6])
			TEMP_vcf_ind[,8] = gsub("CIPOS=","",TEMP_vcf_ind[,8])
			TEMP_vcf_ind[,9] = gsub("CIEND=","",TEMP_vcf_ind[,9])	
		}
		
		if(SVcaller=="Lumpy"){
			tmp_vcf_filter = subset(tmp_vcf,tmp_vcf[,5]=="<DEL>" | tmp_vcf[,5]=="<INS>" | tmp_vcf[,5]=="<DUP>")
			tmp_vcf_1 = tmp_vcf_filter[,c("CHROM","POS","REF","ALT")]
			tmp_vcf_2 = do.call("rbind",strsplit(tmp_vcf_filter[,8],";"))[,c(4,1,3,5,6)]
			tmp_vcf_3 = do.call("rbind",strsplit(tmp_vcf_filter[,10],":"))[,1:4]
			tmp_vcf_4 = cbind(tmp_vcf_1,tmp_vcf_2,tmp_vcf_3)[,c(1,2,5,3,4,6:13)]
			TEMP_vcf_ind = cbind(tmp_vcf_4[,1:13])
			colnames(TEMP_vcf_ind) = c("CHROM","POS","END","REF","ALT","SVTYPE","SVLEN","CIPOS","CIEND","GT","SU","PE","SR")
			TEMP_vcf_ind[,3] = gsub("END=","",TEMP_vcf_ind[,3])
			TEMP_vcf_ind[,6] = gsub("SVTYPE=","",TEMP_vcf_ind[,6])
			TEMP_vcf_ind[,7] = gsub("SVLEN=","",TEMP_vcf_ind[,7])
			TEMP_vcf_ind[,8] = gsub("CIPOS=","",TEMP_vcf_ind[,8])
			TEMP_vcf_ind[,9] = gsub("CIEND=","",TEMP_vcf_ind[,9])	
		}
		
		if(SVcaller=="Manta"){
			tmp_vcf_Manta_check = cbind(as.data.frame(do.call("rbind",strsplit(tmp_vcf[,3],":"))[,1]),tmp_vcf)
			tmp_vcf_filter = subset(tmp_vcf_Manta_check,tmp_vcf_Manta_check[,1]=="MantaDEL" | tmp_vcf_Manta_check[,1]=="MantaINS" | tmp_vcf_Manta_check[,1]=="MantaDUP")[,-1]
			tmp_vcf_1 = tmp_vcf_filter[,c("CHROM","POS","REF","ALT")]
			tmp_vcf_2 = do.call("rbind",strsplit(tmp_vcf_filter[,8],";"))[,c(1,2,3,5,8:10)]
			tmp_vcf_3 = cbind(tmp_vcf_1,tmp_vcf_2)[,c(1,2,5,3,4,6:11)]
			TEMP_vcf_ind = cbind(tmp_vcf_3[,1:8],"NA","./.",tmp_vcf_3[,9:11])
			colnames(TEMP_vcf_ind) = c("CHROM","POS","END","REF","ALT","SVTYPE","SVLEN","CIPOS","CIEND","GT","UPSTREAM_PAIR_COUNT","DOWNSTREAM_PAIR_COUNT","PAIR_COUNT")
			TEMP_vcf_ind[,3] = gsub("END=","",TEMP_vcf_ind[,3])
			TEMP_vcf_ind[,6] = gsub("SVTYPE=","",TEMP_vcf_ind[,6])
			TEMP_vcf_ind[,7] = gsub("SVLEN=","",TEMP_vcf_ind[,7])
			TEMP_vcf_ind[,8] = gsub("CIPOS=","",TEMP_vcf_ind[,8])
			TEMP_vcf_ind[,11] = gsub("UPSTREAM_PAIR_COUNT=","",TEMP_vcf_ind[,11])
			TEMP_vcf_ind[,12] = gsub("DOWNSTREAM_PAIR_COUNT=","",TEMP_vcf_ind[,12])
			TEMP_vcf_ind[,13] = gsub("PAIR_COUNT=","",TEMP_vcf_ind[,13])
			
		}
		
		if(SVcaller=="SoftSearch"){
			tmp_vcf_SoftSearch_check = cbind(as.data.frame(do.call("rbind",strsplit(tmp_vcf[,8],";"))[,2]),tmp_vcf)
			tmp_vcf_filter = subset(tmp_vcf_SoftSearch_check,tmp_vcf_SoftSearch_check[,1]=="EVENT=DEL" | tmp_vcf_SoftSearch_check[,1]=="EVENT=NOV_INS" | tmp_vcf_SoftSearch_check[,1]=="EVENT=TDUP")[,-1]
			tmp_vcf_1 = tmp_vcf_filter[,c("CHROM","POS","REF","ALT")]
			tmp_vcf_2 = do.call("rbind",strsplit(tmp_vcf_filter[,8],";"))[,c(3,2,4)]
			tmp_vcf_3 = do.call("rbind",strsplit(tmp_vcf_filter[,10],":"))[,c(1:11)]
			tmp_vcf_4 = cbind(tmp_vcf_1,tmp_vcf_2,tmp_vcf_3)[,c(1,2,5,3,4,6:18)]
			TEMP_vcf_ind = cbind(tmp_vcf_4[,1:7],"NA","NA",tmp_vcf_4[8:18])
			colnames(TEMP_vcf_ind) = c("CHROM","POS","END","REF","ALT","SVTYPE","SVLEN","CIPOS","CIEND","GT","CTX","DEL","INS","INV","NOV_INS","TDUP","lSC","nSC","uRP","distl_levD")
			TEMP_vcf_ind[,1] = gsub("chr","",TEMP_vcf_ind[,1])
			TEMP_vcf_ind[,3] = gsub("END=","",TEMP_vcf_ind[,3])
			TEMP_vcf_ind[,6] = gsub("EVENT=DEL","DEL",TEMP_vcf_ind[,6])
			TEMP_vcf_ind[,6] = gsub("EVENT=NOV_INS","INS",TEMP_vcf_ind[,6])
			TEMP_vcf_ind[,6] = gsub("EVENT=TDUP","DUP",TEMP_vcf_ind[,6])
			TEMP_vcf_ind[,7] = gsub("ISIZE=","",TEMP_vcf_ind[,7])
		}
		
		if(SVcaller=="SpeedSeq"){
			tmp_vcf_filter = subset(tmp_vcf,tmp_vcf[,5]=="<DEL>" | tmp_vcf[,5]=="<INS>" | tmp_vcf[,5]=="<DUP>")
			tmp_vcf_1 = tmp_vcf_filter[,c("CHROM","POS","REF","ALT")]
			tmp_vcf_2 = do.call("rbind",strsplit(tmp_vcf_filter[,8],";"))[,c(4,1,3,5,6)]
			tmp_vcf_3 = do.call("rbind",strsplit(tmp_vcf_filter[,10],":"))[,1:4]
			tmp_vcf_4 = cbind(tmp_vcf_1,tmp_vcf_2,tmp_vcf_3)[,c(1,2,5,3,4,6:13)]
			TEMP_vcf_ind = cbind(tmp_vcf_4[,1:13])
			colnames(TEMP_vcf_ind) = c("CHROM","POS","END","REF","ALT","SVTYPE","SVLEN","CIPOS","CIEND","GT","SU","PE","SR")
			TEMP_vcf_ind[,3] = gsub("END=","",TEMP_vcf_ind[,3])
			TEMP_vcf_ind[,6] = gsub("SVTYPE=","",TEMP_vcf_ind[,6])
			TEMP_vcf_ind[,7] = gsub("SVLEN=","",TEMP_vcf_ind[,7])
			TEMP_vcf_ind[,8] = gsub("CIPOS=","",TEMP_vcf_ind[,8])
			TEMP_vcf_ind[,9] = gsub("CIEND=","",TEMP_vcf_ind[,9])	
		}
		
		if(SVcaller=="Wham"){
			tmp_vcf_filter = subset(tmp_vcf,tmp_vcf[,5]=="<DEL>" | tmp_vcf[,5]=="<INS>" | tmp_vcf[,5]=="<DUP>")
			tmp_vcf_1 = tmp_vcf_filter[,c("CHROM","POS","REF","ALT")]
			tmp_vcf_2 = do.call("rbind",strsplit(tmp_vcf_filter[,8],";"))[,c(8,14,13,3,2)]
			tmp_vcf_3 = do.call("rbind",strsplit(tmp_vcf_filter[,10],":"))[,1:3]
			tmp_vcf_4 = cbind(tmp_vcf_1,tmp_vcf_2,tmp_vcf_3)[,c(1,2,5,3,4,6:12)]
			TEMP_vcf_ind = cbind(tmp_vcf_4[,1:12])
			colnames(TEMP_vcf_ind) = c("CHROM","POS","END","REF","ALT","SVTYPE","SVLEN","CIPOS","CIEND","GT","DP","SP")
			TEMP_vcf_ind[,3] = gsub("END=","",TEMP_vcf_ind[,3])
			TEMP_vcf_ind[,6] = gsub("SVTYPE=","",TEMP_vcf_ind[,6])
			TEMP_vcf_ind[,7] = gsub("SVLEN=","",TEMP_vcf_ind[,7])
			TEMP_vcf_ind[,8] = gsub("CIPOS=","",TEMP_vcf_ind[,8])
			TEMP_vcf_ind[,9] = gsub("CIEND=","",TEMP_vcf_ind[,9])	
		}
		
		#assign(paste("tmp_vcf_",Ind,sep=""),TEMP_vcf_ind)
		system(paste("rm ",SVcaller,"/tmp_vcf*",sep=""))
		TEMP_vcf_ind[,"POS"] = as.numeric(TEMP_vcf_ind[,"POS"])
		TEMP_vcf_ind[,"END"] = as.numeric(TEMP_vcf_ind[,"END"])
		write.table(TEMP_vcf_ind,file=paste(SVcaller,"/4_MERGED_TXT/",SVcaller,"_",Ind,".txt",sep=""),row.names=F,col.names=T,quote=F)
		cat(paste(SVcaller," : ",Ind,"'s SV-TXT file has been generated","\n",sep=""))
	}else{
		cat(paste(SVcaller," : ",Ind,"'s SV-TXT file has not been generated","\n",sep=""))
	}
}
