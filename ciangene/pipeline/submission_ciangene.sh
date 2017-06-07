release=March2017

repo=/SAN/vyplab/UCLex/scripts/DNASeq_pipeline/ciangene/
export repo
pipeline=${repo}/pipeline/pipeline_cian_gene.sh


script=cluster/submission/cian.sh

rootODir=/SAN/vyplab/UCLex/mainset_${release}/cian/
#rootODir=/cluster/scratch3/vyp-scratch2/cian

step1=yes
step2=no
step3=no
step4=no
step5=no

mkdir -p $rootODir

sh ${pipeline} --step1 ${step1} --step2 ${step2} --step3 ${step3} --step4 ${step4} --step5 ${step5} --rootODir ${rootODir} --release $release
