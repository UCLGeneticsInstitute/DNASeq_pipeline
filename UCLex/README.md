
# Building new UCLex release

The current working directory is here:
```
/SAN/vyplab/UCLex/
```

The main script:
```
msample_calling.sh
```
You can run the individual steps:
```
release=July2016
```
## Step 1: combine gVCF files

Combine the gVCFs:
```
bash ./scripts/DNASeq_pipeline/UCLex/msample_calling.sh --genotype yes --gVCFlist gvcf_list.mainset_${release} --currentUCLex ${release}
```
This creates:
```
 mainset_September2015_chr9.vcf.gz
 mainset_September2015_chr9.vcf.gz.tbi
 ```

Recalibrate the VCFs:
```
bash ./scripts/DNASeq_pipeline/UCLex/msample_calling.sh --recal yes --gVCFlist gvcf_list.mainset_${release} --currentUCLex ${release}
```
Annotate the variants.
The input of which is:
```
 mainset_September2015_chr9_db
 mainset_September2015_chr9_db.log
```
```
bash ./scripts/DNASeq_pipeline/UCLex/msample_calling.sh --annovar yes --gVCFlist gvcf_list.mainset_${release} --currentUCLex ${release}
```
Convert to R datasets:
```
bash ./scripts/DNASeq_pipeline/UCLex/msample_calling.sh --convertToR yes --gVCFlist gvcf_list.mainset_${release} --currentUCLex ${release}
```

 GATK filter and VSQR:
 ```
 mainset_September2015_chr9_filtered.vcf.idx
 mainset_September2015_chr9_indels_filtered.vcf.gz
 mainset_September2015_chr9_indels_filtered.vcf.gz.tbi
 mainset_September2015_chr9_indels.vcf.gz.tbi
 
 mainset_September2015_chr9_recal_plots_snps.R
 mainset_September2015_chr9_recal_plots_snps.R.pdf
 
 mainset_September2015_chr9_SNPs_combrec.recal
 mainset_September2015_chr9_SNPs_combrec.recal.idx
 mainset_September2015_chr9_SNPs_combtranch
 mainset_September2015_chr9_SNPs_filtered.vcf.gz.tbi
 mainset_September2015_chr9_SNPs.vcf.gz.tbi
```



## Step 2: msample_calling.sh: Massive joint calling GenotypeGVCFs

The following are all steps in ``` msample_calling.sh ```:

### GATK GenotypeGVCFs: Joint calling

Input:
```
    if [[ "$format" == "v1" ]]; then gVCF=${path}/chr${chr}/${id}.gvcf.gz; fi
    if [[ "$format" == "v2" ]]; then gVCF=${path}/${id}-chr${chr}.gvcf.gz; fi
    if [[ "$format" == "v3" ]]; then gVCF=${path}/chr${chr}.gvcf.gz; fi
```

Output:
```
${output}_chr${chr}.vcf.gz
```

Code:
```
$java -Xmx${memoSmall}g -jar $GATK -R $fasta \\
-T GenotypeGVCFs \\
-L $chr -L $target --interval_set_rule INTERSECTION --interval_padding 100  \\
--annotation InbreedingCoeff --annotation QualByDepth --annotation HaplotypeScore --annotation MappingQualityRankSumTest --annotation ReadPosRankSumTest --annotation FisherStrand \\
--dbsnp ${bundle}/dbsnp_137.b37.vcf \\" >> cluster/submission/subscript_chr${chr}.sh
while read path id format; do
    if [[ "$format" == "v1" ]]; then gVCF=${path}/chr${chr}/${id}.gvcf.gz; fi
    if [[ "$format" == "v2" ]]; then gVCF=${path}/${id}-chr${chr}.gvcf.gz; fi
    if [[ "$format" == "v3" ]]; then gVCF=${path}/chr${chr}.gvcf.gz; fi
    echo "Including $gVCF"
    if [ ! -s $gVCF ]; then stop "Cannot find $gVCF"; fi
    if [ ! -s $gVCF.tbi ]; then stop "Cannot find $gVCF.tbi"; fi
    echo "   --variant $gVCF \\" >> $script
done < <(tail -n +2 $gVCFlist)
echo "   -o ${output}_chr${chr}.vcf.gz" >> $script
```
    
### Variant recalibration and filtering

Input:
```
${output}_chr${chr}.vcf.gz 
```
Output:
```
${output}_chr${chr}_filtered.vcf
```
    
#### GATK SelectVariants: extract the indels

Output:
```
-o ${output}_chr${chr}_indels.vcf.gz
```
Code:
```
$java  -Djava.io.tmpdir=${tmpDir} -Xmx${memoSmall}g -jar ${GATK} \
     -T SelectVariants \
     -R $fasta \
     -V ${output}_chr${chr}.vcf.gz \
     -selectType INDEL \
     -selectType MIXED \
     -o ${output}_chr${chr}_indels.vcf.gz
```

#### GATK VariantFiltration: apply the filters for the indels
Code:
```
$java -Djava.io.tmpdir=${tmpDir} -Xmx${memoSmall}g -jar ${GATK} \
    -T VariantFiltration \
    -R $fasta \
    -V ${output}_chr${chr}_indels.vcf.gz \
    --filterExpression \"QD < 2.0 || FS > 50.0 || ReadPosRankSum < -20.0\" \
    --filterName \"FAIL\" \
    -o ${output}_chr${chr}_indels_filtered.vcf.gz
```


#### GATK SelectVariants: extract the SNPs

Output:
```
${output}_chr${chr}_SNPs.vcf.gz
``` 
Code:
```
$java  -Djava.io.tmpdir=${tmpDir} -Xmx${memoSmall}g -jar ${GATK} \
     -T SelectVariants \
     -R $fasta \
     -V ${output}_chr${chr}.vcf.gz \
     -selectType SNP \
     -o ${output}_chr${chr}_SNPs.vcf.gz
```

#### GATK VariantRecalibrator: SNPs
Input:
```
${output}_chr${chr}_SNPs.vcf.gz 
```
Output:
```
-recalFile ${output}_chr${chr}_SNPs_combrec.recal \
-tranchesFile ${output}_chr${chr}_SNPs_combtranch \
-rscriptFile  ${output}_chr${chr}_recal_plots_snps.R
```
Code:
```
$java -Djava.io.tmpdir=${tmpDir} -Xmx${memoSmall}g -jar ${GATK} -R $fasta 
    -T VariantRecalibrator
    -L $chr --input ${output}_chr${chr}_SNPs.vcf.gz --maxGaussians ${maxGauLoc} --mode SNP \
    -resource:hapmap,VCF,known=false,training=true,truth=true,prior=15.0 ${bundle}/hapmap_3.3.b37.vcf  \
    -resource:omni,VCF,known=false,training=true,truth=false,prior=12.0 ${bundle}/1000G_omni2.5.b37.vcf \
    -resource:dbsnp,VCF,known=true,training=false,truth=false,prior=8.0 ${bundle}/dbsnp_137.b37.vcf \
    -an QD -an FS -an ReadPosRankSum -an InbreedingCoeff \
    -tranche 100.0 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 \
    --minNumBadVariants ${numBad} \
    -recalFile ${output}_chr${chr}_SNPs_combrec.recal \
    -tranchesFile ${output}_chr${chr}_SNPs_combtranch \
    -rscriptFile  ${output}_chr${chr}_recal_plots_snps.R
```

Make the R plots
```
${Rscript} ${output}_chr${chr}_recal_plots_snps.R
```

#### GATK ApplyRecalibration: recalibration of SNPs

Input:
```
${output}_chr${chr}_SNPs.vcf.gz
```
Output:
```
${output}_chr${chr}_SNPs_filtered.vcf.gz
```

Code:
```
$java -Xmx${memoSmall}g -jar ${GATK} -R $fasta
-T ApplyRecalibration 
-o ${output}_chr${chr}_SNPs_filtered.vcf.gz \
--ts_filter_level 99.5 \
--recal_file ${output}_chr${chr}_SNPs_combrec.recal --tranches_file ${output}_chr${chr}_SNPs_combtranch --mode SNP \
--input ${output}_chr${chr}_SNPs.vcf.gz
```


#### GATK CombineVariants: Now we merge SNPs and indels

Input:
```
--variant:SNPs ${output}_chr${chr}_SNPs_filtered.vcf.gz \
--variant:indels ${output}_chr${chr}_indels_filtered.vcf.gz \
```
Output:
```
${output}_chr${chr}_filtered.vcf
```

Code:
```
$java -Djava.io.tmpdir=${tmpDir} -Xmx${memoSmall}g -jar ${GATK} \
       -T CombineVariants --assumeIdenticalSamples \
       -R $fasta \
       --variant:SNPs ${output}_chr${chr}_SNPs_filtered.vcf.gz \
       --variant:indels ${output}_chr${chr}_indels_filtered.vcf.gz \
       -genotypeMergeOptions PRIORITIZE  \
       -priority SNPs,indels \
       -o ${output}_chr${chr}_filtered.vcf
```

Remove intermediate files:
```
rm ${output}_chr${chr}_indels.vcf.gz ${output}_chr${chr}_SNPs.vcf.gz ${output}_chr${chr}_SNPs_filtered.vcf.gz
```

## Custom filtering

Perl script to print out as missing, dodgy looking variants.
```
perl ${baseFolder}/GATK_v2/custom_filtering.pl ${output}_chr${chr}_filtered.vcf ${output}_chr${chr}_recal_filtered2.vcf ${GQ}
```

Input:
```
${output}_chr${chr}_filtered.vcf 
```
Output:
```
${output}_chr${chr}_recal_filtered2.vcf
```
If good is false will print as missing:
```
 if ($geno eq "0\/0") {
````````if ($GQ <= $locth) {$good = 0;}
````    } else {
````````my @PLs = split(',', $PL);
````````my $PL0 = $PLs[ 0 ];
````````if ($PL0 <= $locth) {$good = 0;}
````    }

````    #if ($geno eq "1\/1") {$locth = 9;}  ## for ALT/ALT calls, a 15 QUAL will do in any case
````    #if ($geno eq "2\/2") {$locth = 9;}  ## for ALT/ALT calls, a 15 QUAL will do in any case
````    #if ($geno eq "3\/3") {$locth = 9;}  ## for ALT/ALT calls, a 15 QUAL will do in any case

````    if ($ezSNP) {  ##if nice clean di-allelic SNP we add an extra filter
````````if ($geno eq "0\/1") {
````````    my @counts = split(',', $AD);
````````    my $total = $counts[0] + $counts[1];
````````    if ($counts[1] > $counts[0] + 1) {$locth = 9;}  ##if the alternate count is higher then we can be relax, because the alternative is probably

````````    if ($total >= 10) {
````````````my $ratio = $counts[1]/$total;
````````````my $limit = 0.18;
````````````if ($ratio <= $limit) {$good = 0;}
````````    }
````````}
}
```

## Annotation with Annovar

```
${output}_chr${chr}_recal_filtered2.vcf
```

```
${output}_chr${chr}_exome_table.csv 
```

Code:
```
cut -f1-8 ${output}_chr${chr}_filtered.vcf > ${output}_chr${chr}_for_annovar.vcf
/cluster/project8/vyp/vincent/Software_heavy/annovar_Feb2013/convert2annovar.pl --allallele -format vcf4 --includeinfo ${output}_chr${chr}_for_annovar.vcf > ${output}_chr${chr}_db
```

```
/cluster/project8/vyp/vincent/Software_heavy/annovar_Feb2013/summarize_annovar_VP.pl -ver1000g 1000g2012apr -verdbsnp 137 -veresp 6500si -alltranscript -buildver hg19 --genetype ensgene --remove ${output}_chr${chr}_db /cluster/project8/vyp/vincent/Software_heavy/annovar_Feb2013/humandb_hg19/
```

Custom filtering.

Input:
```
${output}_chr${chr}_filtered.vcf
```
Output:
```
${output}_chr${chr}_recal_filtered2.vcf 
```
Code:
```
perl ${baseFolder}/GATK_v2/custom_filtering.pl ${output}_chr${chr}_filtered.vcf ${output}_chr${chr}_recal_filtered2.vcf ${GQ}
```

Combine ANNOVAR and VCF.

Input:
```
${output}_chr${chr}_recal_filtered2.vcf
${output}_chr${chr}_db.exome_summary.csv
```
Output:
```
${output}_chr${chr}_exome_table.csv
```
Code:
```
python ${baseFolder}/GATK_v2/annovar_vcf_combine_VP.py ${output}_chr${chr}_recal_filtered2.vcf ${output}_chr${chr}_db.exome_summary.csv ${output}_chr${chr}_exome_table.csv
```

Extract genotypes:
Input:
```
${output}_chr${chr}_exome_table.csv
```
Output:
```
my $callFile = $output."_calls.tab";
my $depthFile = $output."_depth.tab";
my $colnameFile = $output."_colname.tab";
my $rownameFile = $output."_rowname.tab";
my $annotationFile = $output."_annotations.csv";
```
Code:
```
perl ${baseFolder}/msample_calling/make_matrix_calls.pl ${output}_chr${chr}_exome_table.csv ${output} $chr
```

# convertToR
    
Output:
```
${output}_snpStats/chr${chr}_snpStats.RData
```
Code:
```
R CMD BATCH --no-save --no-restore --chromosome=${chr} --root=${output} ${baseFolder}/msample_calling/convert_to_R.R cluster/R/convert_to_R_chr${chr}.out
```


# crunch_controls.pl: finalCrunch

What's a crunch?

```
perl ${baseFolder}/GATK_v2/crunch_controls.pl ${output}_chr${chr}_exome_table.csv $keyWords $casekeyWords ${output}_chr${chr}_exome_crunched.csv data/sampleList_exome.tab none no  ##include all samples
```

Input:
```
keyWords=data/controlKeywords.tab
casekeyWords=data/caseKeywords.tab
${output}_chr${chr}_exome_table.csv 
```
Output:
```
 ${output}_chr${chr}_exome_crunched.csv 
```

Include all samples:
```
$inputTable = $ARGV[0];
$controlKeywords = $ARGV[1];
$caseKeywords = $ARGV[2];
$outputTable = $ARGV[3];
$sampleList = $ARGV[4];
$excludedControlList = $ARGV[5];
$chooseExternalControlsRandom = $ARGV[6];
crunch_controls.pl ${output}_chr${chr}_exome_table.csv data/controlKeywords.tab data/caseKeywords.tab ${output}_chr${chr}_exome_crunched.csv data/sampleList_exome.tab none no
```

## crunch_data.R: Extract cohort specific information


Output:
```
save(list = c('combined.matrix.depth', 'combined.snpStats', 'combined.annotations'), file = 'data/poly_in_cases_snpStats.RData')
```

Code:
```
crunch_data.R
```

## Finalize the user friendly folders

process_multiVCF.R

# Samples exclusion

sampleExclusion=/cluster/project8/vyp/exome_sequencing_multisamples/mainset/support/exclusion_lists/control_exclusion.tab

# Merge of VCFs

```
head -300 ${output}_chr1 | awk '{ if (\$1 ~ /^#/) print}' > ${output}.vcf
	if [ ! -e ${output}_chr${chr} ]; then
	    echo "Missing file ${output}_chr${chr}"
	else
	    echo "awk '{ if (\$1 !~ /^#/) print}' ${output}_chr${chr} >> ${output}.vcf" >> $mainScript
	fi
    done

fi
```



