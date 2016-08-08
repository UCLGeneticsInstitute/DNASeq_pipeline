#$ -S /bin/bash
#$ -o /dev/null
#$ -e /dev/null
#$ -cwd
#$ -V
#$ -l tmem=14.9G,h_vmem=14.9G
#$ -l h_rt=240:0:0
#$ -t 1-24
set -u
set -x
scriptname=annotate_${UCLEX_RELEASE}

mkdir -p ${scriptname}.${UCLEX_RELEASE}.qsub.out ${scriptname}.${UCLEX_RELEASE}.qsub.err
exec >${scriptname}.${UCLEX_RELEASE}.qsub.out/${scriptname}_${SGE_TASK_ID}_${JOB_ID}.out 2>${scriptname}.${UCLEX_RELEASE}.qsub.err/${scriptname}_${SGE_TASK_ID}_${JOB_ID}.err
args=( header 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y )
f=${args[$SGE_TASK_ID]}
chrCode=$f

#format=pileup
format=vcf
echo SCRATCH2: $SCRATCH2
Software="$SCRATCH2/Software/DNASeq_pipeline/"
mkdir -p $UCLEX_RELEASE
cd $UCLEX_RELEASE
#input=${UCLEX}/${UCLEX_RELEASE}_chr${chrCode}_recal_filtered2.vcf
#input=${UCLEX}/${UCLEX_RELEASE}_chr${chrCode}_for_annovar.vcf
input=${UCLEX}/${UCLEX_RELEASE}_chr${chrCode}_filtered.vcf


function freq() {
    #input=${uclex_basedir}/${UCLEX_RELEASE}_chr${chrCode}.vcf.gz
    #will set to missing individuals with depth < 20
    #if variants are missing in all individuals they are skipped
    # these freq need to be in unrelated individuals
    cat ${input} | python $SCRATCH2/Software/vcf_scripts/vcftools.py --DP 1 --freq > chr${chrCode}.freq
}

#
function VEP() {
    reference=1kg
    #will set to missing individuals with depth < 20
    #if variants are missing in all individuals they are skipped
    #DIR=$SCRATCH2/Software/vcf_scripts/
    #zcat ${input} | python ${Software}/vcf_scripts/vcftools.py --DP 20 --${format} > chr${chrCode}.vcf
    #input=chr${chrCode}.${format}
    DIR="${Software}/annotation"
    cat $input | python ${DIR}/multiallele_to_single_gvcf.py  > chr${chrCode}-single.vcf
    ####CONFIGURE SOFTWARE SHORTCUTS AND PATHS
    ensembl=/cluster/project8/vyp/AdamLevine/software/ensembl/
    #VEP=${ensembl}/src/ensembl-tools/scripts/variant_effect_predictor/variant_effect_predictor.pl
    #dir_cache=${ensembl}/cache/
    PERL5LIB=${PERL5LIB}:${ensembl}/src/bioperl-1.6.1
    PERL5LIB=${PERL5LIB}:${ensembl}/src/ensembl/modules
    PERL5LIB=${PERL5LIB}:${ensembl}/src/ensembl-compara/modules
    PERL5LIB=${PERL5LIB}:${ensembl}/src/ensembl-variation/modules
    PERL5LIB=${PERL5LIB}:${ensembl}/src/ensembl-funcgen/modules
    PERL5LIB=${PERL5LIB}:${ensembl}/Plugins
    export PERL5LIB
    export PATH=$PATH:/cluster/project8/vyp/vincent/Software/tabix-0.2.5/
    /share/apps/perl/bin/perl /cluster/project8/vyp/Software/ensembl-tools-release-82/scripts/variant_effect_predictor/variant_effect_predictor.pl --port 3337 --verbose --ASSEMBLY GRCh37 --fasta /cluster/scratch3/vyp-scratch2//reference_datasets/human_reference_sequence//human_g1k_v37.fasta --cache --dir_cache /cluster/project8/vyp/AdamLevine/software/ensembl//cache/ --sift b --polyphen b --symbol --canonical --check_existing --check_alleles --no_progress --fork 4 --maf_esp --gmaf --maf_1kg \
    --custom /cluster/project9/IBDAJE/VEP_custom_annotations/1kg/CADD/chr${chr}.vcf.gz,CADD,vcf,exact \
    --custom /cluster/project9/IBDAJE/VEP_custom_annotations/1kg/ExAC/0.3/chr${chr}_AFR.vcf.gz,EXAC_AFR,vcf,exact \
    --custom /cluster/project9/IBDAJE/VEP_custom_annotations/1kg/ExAC/0.3/chr${chr}_AMR.vcf.gz,EXAC_AMR,vcf,exact \
    --custom /cluster/project9/IBDAJE/VEP_custom_annotations/1kg/ExAC/0.3/chr${chr}_Adj.vcf.gz,EXAC_Adj,vcf,exact \
    --custom /cluster/project9/IBDAJE/VEP_custom_annotations/1kg/ExAC/0.3/chr${chr}_EAS.vcf.gz,EXAC_EAS,vcf,exact \
    --custom /cluster/project9/IBDAJE/VEP_custom_annotations/1kg/ExAC/0.3/chr${chr}_FIN.vcf.gz,EXAC_FIN,vcf,exact \
    --custom /cluster/project9/IBDAJE/VEP_custom_annotations/1kg/ExAC/0.3/chr${chr}_NFE.vcf.gz,EXAC_NFE,vcf,exact \
    --custom /cluster/project9/IBDAJE/VEP_custom_annotations/1kg/ExAC/0.3/chr${chr}_OTH.vcf.gz,EXAC_OTH,vcf,exact \
    --custom /cluster/project9/IBDAJE/VEP_custom_annotations/1kg/ExAC/0.3/chr${chr}_SAS.vcf.gz,EXAC_SAS,vcf,exact \
    --custom /cluster/project9/IBDAJE/VEP_custom_annotations/1kg/1kg/chr${chr}_EUR.vcf.gz,1KG_EUR,vcf,exact \
    --custom /cluster/project9/IBDAJE/VEP_custom_annotations/1kg/1kg/chr${chr}_AFR.vcf.gz,1KG_AFR,vcf,exact \
    --custom /cluster/project9/IBDAJE/VEP_custom_annotations/1kg/1kg/chr${chr}_AMR.vcf.gz,1KG_AMR,vcf,exact \
    --custom /cluster/project9/IBDAJE/VEP_custom_annotations/1kg/1kg/chr${chr}_ASN.vcf.gz,1KG_ASN,vcf,exact \
    --custom /cluster/project9/IBDAJE/VEP_custom_annotations/1kg/esp/chr${chr}_EA.vcf.gz,ESP_EA,vcf,exact \
    --custom /cluster/project9/IBDAJE/VEP_custom_annotations/1kg/esp/chr${chr}_AA.vcf.gz,ESP_AA,vcf,exact \
    --custom /cluster/scratch3/vyp-scratch2/reference_datasets/Kaviar/Kaviar-160204-Public/hg19/VEP_annotation.vcf.gz,Kaviar,vcf,exact
    --plugin Condel,/cluster/project8/vyp/AdamLevine/software/ensembl//Plugins/config/Condel/config,b \
    --plugin Carol \
    --plugin CADD,/cluster/project9/IBDAJE/VEP_custom_annotations/1kg/CADD/chr${chr}.vcf.gz \
    --plugin GO \
    --plugin ExAC,/cluster/project9/IBDAJE/VEP_custom_annotations/1kg/ExAC/0.3/chr${chr}.vcf.gz \
    --input_file chr${chrCode}-single.vcf \
    --force_overwrite \
    --hgvs \
    --plugin HGVSshift \
    --plugin LoFtool,scores_file.txt \
    --plugin SameCodon \
    --json --output_file VEP_${chr}.json
}

#
function prepare_for_uclex() {
   out="$SCRATCH2/vincent/uclex_browser_data/${UCLEX_RELEASE}_chr${chrCode}.vcf"
   python $Software/uclex_browser/uclex_to_vcf.py --infile VEP_${chrCode}.vcfout >  $out
   rm -f ${out}.gz*
   bgzip ${out}
   tabix -f -p vcf ${out}.gz
}

function to_json() {
    bgziptabix VEP_${chrCode}.vcfout
    python /home/rmhanpo/uclex_browser/load_vcfs/vcf_to_json.py --infile VEP_${chrCode}.vcfout | gzip -c  > VEP_${chrCode}.json.gz
    scp VEP_${chrCode}.json.gz phenotips:/slms/UGI/vm_exports/vyp/phenotips/uclex_files/$UCLEX_RELEASE/
}


#freq
VEP
#annotate
#prepare_for_uclex
#to_json


