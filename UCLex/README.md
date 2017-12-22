
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
release=November2017
```

## Preliminary steps for new exomes: align reads, make gvcfs, combine gvcfs.

First make a file with 3 white-space separated columns called 'code', 'f1' and 'f2' which contain the sample names and fastq file paths. Prefix sample names with the projectID e.g. BGI_Nov2017_116samples_IBDAJ. Sitting in the same directory, create variables called 'projectID' and 'samples' and set them to the projectID and samples file respectively. Then run the following commands to make and submit the alignment job scripts:
```
bash /SAN/vyplab/UCLex/scripts/DNASeq_pipeline/WGS/WGS_pipeline.sh --mode align --supportFrame $samples --reference 1kg --aligner-tparam 320 --inputFormat STDFQ --projectID $projectID --vmem 1
qsub $projectID/align/scripts/align.sh
```
Bam files are written to /SAN/vyplab/UCLex_raw. If jobs fail, then re-run the bash command to make still outstanding jobs. Next, to make and submit the job scripts for making the gvcfs, run:
```
bash /SAN/vyplab/UCLex/scripts/DNASeq_pipeline/WGS/WGS_pipeline.sh --mode gvcf --supportFrame $samples --reference 1kg --inputFormat STDFQ --projectID  $projectID --vmem 14
qsub $projectID/gvcf/scripts/gvcf.sh
```
The gvcfs are written to $projectID/gvcf/data. Finally, to make and submit the job scripts for combining the gvcf, run:
```
bash /SAN/vyplab/UCLex/scripts/DNASeq_pipeline/WGS/WGS_pipeline.sh --mode CombineGVCFs --supportFrame $samples --reference 1kg --inputFormat STDFQ --projectID $projectID --batchName $projectID
qsub $projectID/CombineGVCFs/scripts/CombineGVCFs.sh
```
The combined gvcfs are written to $UCLEX_DIR/combinedGVCFs i.e. /SAN/vyplab/UCLex/combinedGVCFs.

## Step 1: combine gVCF files

First make a new gvcf_list.mainset file in $UCLEX_DIR that contains all batches for the new build i.e. batches in the previous build plus the new ones. Use ```./scripts/DNASeq_pipeline/UCLex/msample_calling.sh``` under different modes to create job scripts:
```
bash ./scripts/DNASeq_pipeline/UCLex/msample_calling.sh  --gVCFlist gvcf_list.mainset_${release} --currentUCLex ${release} --mode <MODE>
```
then submit them with ```qsub $UCLEX_DIR/mainset_${release}/scripts/calling.sh```. First combine the gVCF files using mode "genotype".

## Step 2: recalibrate the SNPs and InDels

Extract the SNPs: mode "extract_snps".</br>
Recalibrate the SNPs: mode "recal_snps".</br>
Extract the InDels: mode "extract_indels".</br>
Recalibrate the InDels: mode "recal_indels".</br>
Merge the recalibrated SNPs and InDels: mode "recal_merge".</br>

## Step 3: ANNOVAR

Run ANNOVAR: mode "annovar".</br>
Convert to R datasets: mode "convertToR".</br>

### Create giant VCF file

This will concatenate all chromosome VCF files and store it in the bgt folder:
```
cd ${UCLEX}
bash concat.sh
```

### BGT

[BGT](https://academic.oup.com/bioinformatics/article/32/4/590/1743991) is a efficient tool for querying variant calls.
We create the BGT file on the giant VCF:
```
cd ${UCLEX}
bash bgt.sh
```

### Split multi-allelelic variants across lines

This is required for CADD and VEP.
```
cd ${UCLEX}
bash multi2single.sh
```

### CADD

```
cd ${UCLEX}
bash cadd.sh
```

## Variant Effect Predictor (VEP)

VEP annotation now runs on beck server maintained by @IsmailM.
Files need to prepared for VEP and then transferred to the server.
The files are:
* for_VEP.vcf (created by multi2single.sh)
* CADD.vcf.gz (created by cadd.sh)
VEP outputs JSON that can then be processed by the script written by @mframpton.

Make input for VEP: mode "VEP_input".</br>
Get CADD scores: mode "CADD".</br>
Run script to make for_vep.vcf.</br>
Feed for_vep.vcf to VEP and get back json output.</br>
Run split_json_vep_by_chr.py --uclex_build November2017 to split the json output into per-chromosome csvs.</br>

## End of November2017 updates

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



