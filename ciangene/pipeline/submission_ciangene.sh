#repo=/cluster/project8/vyp/exome_sequencing_multisamples/ciangene/
repo=/SAN/vyplab/UCLex/scripts/DNASeq_pipeline/ciangene/

export repo

pipeline=${repo}/pipeline/pipeline_cian_gene.sh


script=cluster/submission/cian.sh
rootODir=/SAN/vyplab/UCLex/mainset_July2016/cian/
#rootODir=/cluster/scratch3/vyp-scratch2/cian

step1=yes
step2=no
step3=no
step4=no 

mkdir -p $rootODir

sh ${pipeline} --step1 ${step1} --step2 ${step2} --step3 ${step3} --step4 ${step4} --rootODir ${rootODir} --release July2016
