source('/SAN/vyplab/UCLex/scripts/DNASeq_pipeline/ciangene/Gene_Based_Tests/SKAT/summarise.results.function.R')
dbup='/cluster/project8/vyp/cian/scripts/bash/dropbox_uploader.sh'
######################################################################################################################################
######################################################################################################################################
mac<-FALSE

if(mac)
{
genes='/SAN/vyplab/UCLex/support/IoO/Macular_Dystrophy/macular_dystrophy_genes'
summarise(dir='/SAN/vyplab/UCLex/mainset_July2016/cian/HPO/Macular_Dystrophy_cpd_hets',
	genes,Title='Macular_Dystrophy_cpd_hets',
	outputDirectory='/SAN/vyplab/UCLex/mainset_July2016/cian/HPO/Results/Macular_Dystrophy_cpd_hets')
}

######################################################################################################################################
genes='/SAN/vyplab/UCLex/support/Sisodiya/seizures-epilepsy_candidate-gene-list.txt'

oDir='/SAN/vyplab/UCLex/mainset_July2016/cian/HPO/Results/Epilepsy_no_known_cause/'
if(file.exists(oDir))file.remove(oDir)
summarise(dir='/SAN/vyplab/UCLex/mainset_July2016/cian/HPO/Epilepsy/No_known_cause_unrelated/',
	genes,Title='Epilepsy_no_known_cause',
	outputDirectory=oDir)

run<-paste(dbup,'upload',oDir, paste0('PostDoc/Epilepsy/',basename(oDir))) 
system(run)


oDir='/SAN/vyplab/UCLex/mainset_July2016/cian/HPO/Results/Epilepsy_all_cases/'
if(file.exists(oDir))file.remove(oDir)
summarise(dir='/SAN/vyplab/UCLex/mainset_July2016/cian/HPO/Epilepsy/unrelated',
	genes,Title='Epilepsy_all_cases',
	outputDirectory=oDir)
run<-paste(dbup,'upload',oDir, paste0('PostDoc/Epilepsy/',basename(oDir))) 
system(run)


######################################################################################################################################
genes='/SAN/vyplab/UCLex/support/IoO/Macular_Dystrophy/macular_dystrophy_genes'

oDir='/SAN/vyplab/UCLex/mainset_July2016/cian/HPO/Results/Macular_Dystrophy'
if(file.exists(oDir))file.remove(oDir)
summarise(dir='/SAN/vyplab/UCLex/mainset_July2016/cian/HPO/Macular_Dystrophy',
	genes,Title='Macular_Dystrophy',
	outputDirectory=oDir)
run<-paste(dbup,'upload',oDir, paste0('PostDoc/IoO/',basename(oDir))) 
system(run)

oDir='/SAN/vyplab/UCLex/mainset_July2016/cian/HPO/Results/occult_macular_dystrophy'
if(file.exists(oDir))file.remove(oDir)
summarise(dir='/SAN/vyplab/UCLex/mainset_July2016/cian/HPO/occult_macular_dystrophy',
	genes,Title='occult_macular_dystrophy',
	outputDirectory=oDir)
run<-paste(dbup,'upload',oDir, paste0('PostDoc/IoO/',basename(oDir))) 
system(run)

oDir='/SAN/vyplab/UCLex/mainset_July2016/cian/HPO/Results/OFG'
if(file.exists(oDir))file.remove(oDir)
summarise(dir='/SAN/vyplab/UCLex/mainset_July2016/cian/HPO/OFG',
	genes,Title='OFG',
	outputDirectory=oDir)
run<-paste(dbup,'upload',oDir, paste0('PostDoc/IoO/',basename(oDir))) 
system(run)

######################################################################################################################################

oDir='/SAN/vyplab/UCLex/mainset_July2016/cian/HPO/Results/Vulliamy_Bone_Marrow_Failure'
if(file.exists(oDir))file.remove(oDir)
summarise(dir='/SAN/vyplab/UCLex/mainset_July2016/cian/HPO/Vulliamy',
	genes,Title='Vulliamy_Bone_Marrow_Failure',
	outputDirectory=oDir)
run<-paste(dbup,'upload',oDir, paste0('PostDoc/',basename(oDir))) 
system(run)
######################################################################################################################################

genes='/cluster/project8/vyp/cian/data/Support/CandidateGenes/SCD.candidate.genes'

oDir='/SAN/vyplab/UCLex/mainset_July2016/cian/HPO/Results/SCD_Batch1'
if(file.exists(oDir))file.remove(oDir)
summarise(dir='/SAN/vyplab/UCLex/mainset_July2016/cian/HPO/SCD/Batch1',
	genes,Title='SCD_Batch1',
	outputDirectory=oDir)
run<-paste(dbup,'upload',oDir, paste0('PostDoc/SCD/',basename(oDir))) 
system(run)

oDir='/SAN/vyplab/UCLex/mainset_July2016/cian/HPO/Results/SCD_Batch2'
if(file.exists(oDir))file.remove(oDir)
summarise(dir='/SAN/vyplab/UCLex/mainset_July2016/cian/HPO/SCD/Batch2',
	genes,Title='SCD_Batch2',
	outputDirectory=oDir)
run<-paste(dbup,'upload',oDir, paste0('PostDoc/SCD/',basename(oDir))) 
system(run)

oDir='/SAN/vyplab/UCLex/mainset_July2016/cian/HPO/Results/VF'
if(file.exists(oDir))file.remove(oDir)
summarise(dir='/SAN/vyplab/UCLex/mainset_July2016/cian/HPO/SCD/VF',
	genes,Title='VF',
	outputDirectory='/SAN/vyplab/UCLex/mainset_July2016/cian/HPO/Results/VF')
run<-paste(dbup,'upload',oDir, paste0('PostDoc/SCD/',basename(oDir))) 
system(run)