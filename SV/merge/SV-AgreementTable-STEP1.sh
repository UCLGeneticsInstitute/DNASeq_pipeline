
#############################
### STEP1: Copy vcf files ###
#############################
cp /SAN/vyplab/NCMD/b37/SV/Delly/DEL/DEL2.bcf* /SAN/vyplab/NCMD/b37/SV-MERGED/Delly/1_ORIGINAL
cp /SAN/vyplab/NCMD/b37/SV/Delly/INS/INS2.bcf* /SAN/vyplab/NCMD/b37/SV-MERGED/Delly/1_ORIGINAL

cp /SAN/vyplab/NCMD/b37/SV/Lumpy/*.vcf /SAN/vyplab/NCMD/b37/SV-MERGED/Lumpy/1_ORIGINAL

for file in $(ls -d -1 /SAN/vyplab/NCMD/b37/SV/Manta/FR* | perl -pe 's/\n/ /g')
do
 x=$file; z=${x##*/}
 cp /SAN/vyplab/NCMD/b37/SV/Manta/$z/results/variants/candidateSV.vcf.gz /SAN/vyplab/NCMD/b37/SV-MERGED/Manta/1_ORIGINAL/$z.vcf.gz
 cp /SAN/vyplab/NCMD/b37/SV/Manta/$z/results/variants/candidateSV.vcf.gz.tbi /SAN/vyplab/NCMD/b37/SV-MERGED/Manta/1_ORIGINAL/$z.vcf.gz.tbi
 echo $z
done

cp /SAN/vyplab/NCMD/b37/SV/SoftSearch/FR*/*.vcf /SAN/vyplab/NCMD/b37/SV-MERGED/SoftSearch/1_ORIGINAL

cp /SAN/vyplab/NCMD/b37/SV/SpeedSeq/*.vcf.gz* /SAN/vyplab/NCMD/b37/SV-MERGED/SpeedSeq/1_ORIGINAL

cp /SAN/vyplab/NCMD/b37/SV/Wham/WHAMG-CHRS/Chr*/WHAMG.Chr*.vcf /SAN/vyplab/NCMD/b37/SV-MERGED/Wham/1_ORIGINAL

#####################################################################
### STEP2: From vcf to vcf.gz (Lumpy & SoftSearch & Wham) ###
#####################################################################
cd /SAN/vyplab/NCMD/b37/SV-MERGED/Lumpy/1_ORIGINAL
export PERL5LIB=$PERL5LIB:/home/ucbtlcu/bin/0-vcftools_0.1.13/lib/perl5/site_perl/ # in case the vcftools report error #

for file in *.vcf;
do
 x=$file
 y=${x%.vcf} # this removes ".vcf" #
 z=${y##*/} # this removes everything before the final / #
 cat $z.vcf | vcf-sort > ../2_SEPARATED/$z.sorted.vcf
 bgzip ../2_SEPARATED/$z.sorted.vcf 
 tabix -p vcf ../2_SEPARATED/$z.sorted.vcf.gz
 echo $z
done

###########################################################################################
### STEP3-1: Merge IND-separated files together (Lumpy & Manta & SoftSearch & SpeedSeq) ###
###########################################################################################
cd /SAN/vyplab/NCMD/b37/SV-MERGED/Lumpy/2_SEPARATED

# Method1: Using bcftools merge, http://www.htslib.org/doc/bcftools.html#merge #
bcftools merge $(ls -1 *.vcf.gz | perl -pe 's/\n/ /g') -O z -o ../3_MERGED/INDs_merged1.vcf.gz
tabix -p vcf ../3_MERGED/INDs_merged1.vcf.gz

# Method2: Using vcftools vcf-merge, http://vcftools.sourceforge.net/perl_module.html#vcf-merge #
vcf-merge $(ls -1 *.vcf.gz | perl -pe 's/\n/ /g') | bgzip -c > ../3_MERGED/INDs_merged2.vcf.gz
tabix -p vcf ../3_MERGED/INDs_merged2.vcf.gz

# Method3: Using GATK CombineVariants, https://software.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_variantutils_CombineVariants.php #
/share/apps/java/bin/java -Xms2g -jar /home/ucbtlcu/bin/picard.jar CreateSequenceDictionary R=/SAN/vyplab/UKIRDC/reference/human_g1k_v37.fasta O=/SAN/vyplab/UKIRDC/reference/human_g1k_v37.dict
ls -1 *.vcf.gz > input.list
/share/apps/java/bin/java -Xms2g -jar /home/ucbtlcu/bin/GenomeAnalysisTK.jar -T CombineVariants -R /SAN/vyplab/UKIRDC/reference/human_g1k_v37.fasta --variant input.list -o ../3_MERGED/INDs_merged3.vcf -genotypeMergeOptions UNIQUIFY -nt 4
bgzip ../3_MERGED/INDs_merged3.vcf 
tabix -p vcf ../3_MERGED/INDs_merged3.vcf.gz


#########################################################################
### STEP3-2: Merge CHR-separated vcf.gz files together (Delly & Wham) ###
#########################################################################
cd /SAN/vyplab/NCMD/b37/SV-MERGED/Wham/2_SEPARATED

# Method1: Using Shell generate manually # 
bcftools view -h WHAMG.Chr10.sorted.vcf.gz >> ../3_MERGED/CHRs_merged1.vcf # Just watch out to change the file name #
for file in $(ls -1 *.vcf.gz | perl -pe 's/\n/ /g')
do
  bcftools view -H $file >> ../3_MERGED/CHRs_merged1.vcf
  echo $file
done
cat ../3_MERGED/CHRs_merged1.vcf | vcf-sort > ../3_MERGED/CHRs_merged1.sorted.vcf
bgzip ../3_MERGED/CHRs_merged1.sorted.vcf
tabix -p vcf ../3_MERGED/CHRs_merged1.sorted.vcf.gz

# Method2: Using vcftools vcf-concat, http://vcftools.sourceforge.net/perl_module.html#vcf-concat #
vcf-concat $(ls -1 *.vcf.gz | perl -pe 's/\n/ /g') > ../3_MERGED/CHRs_merged2.vcf
cat ../3_MERGED/CHRs_merged2.vcf | vcf-sort > ../3_MERGED/CHRs_merged2.sorted.vcf
bgzip ../3_MERGED/CHRs_merged2.sorted.vcf
tabix -p vcf ../3_MERGED/CHRs_merged2.sorted.vcf.gz


###############################################################################################
### STEP3-3: Separate the WHOLE vcf files into different IND-separated files (Delly & Wham) ###
###############################################################################################
cd /SAN/vyplab/NCMD/b37/SV-MERGED/Wham/3_MERGED

# Method1: Using bcftools view # 
for file in $(bcftools view -h CHRs_merged1.sorted.vcf.gz | tail -1 | awk '{for( i=10; i<=NF; i++ ){printf( "%s ", $i )}}')
do
 bcftools view -s $file CHRs_merged1.sorted.vcf.gz -O z -o ../2_SEPARATED/$file.sorted.vcf.gz
 tabix -p vcf ../2_SEPARATED/$file.sorted.vcf.gz
 echo $file
done

# Method2: Using vcftools vcf-subset #
for file in $(bcftools view -h CHRs_merged1.sorted.vcf.gz | tail -1 | awk '{for( i=10; i<=NF; i++ ){printf( "%s ", $i )}}')
do
 vcf-subset -c $file CHRs_merged1.sorted.vcf.gz | bgzip -c > ../2_SEPARATED/$file.sorted.vcf.gz
 tabix -p vcf ../2_SEPARATED/$file.sorted.vcf.gz
 echo $file
done

#############################################################################################################################################
### STEP4: Using R load IND-separated vcf files from 6 SV-Callers separately, calculate the agreement probability and generate two table. ###
### (1.INDELs from all IND-separated vcf files from 6 SV-Callers; 2.INDELs from all IND-separated vcf files from DELLY)                  ###
#############################################################################################################################################

#source("http://bioconductor.org/biocLite.R") 
#biocLite("VariantAnnotation") #install the package
#library("VariantAnnotation") #load the package
#example("readVcf") #optional, test the function by running example codes
