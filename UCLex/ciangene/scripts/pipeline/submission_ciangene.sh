#repo=/cluster/project8/vyp/exome_sequencing_multisamples/ciangene/
repo=/cluster/project8/vyp/cian/data/UCLex/ciangene/

export repo

pipeline=${repo}/scripts/pipeline/pipeline_cian_gene.sh


script=cluster/submission/cian.sh
# rootODir=/scratch2/vyp-scratch2/ciangene
rootODir=/scratch2/vyp-scratch2/cian/

step1=yes
step2=yes
step3=yes
step4=no 

sh ${pipeline} --step1 ${step1} --step2 ${step2} --step3 ${step3} --step4 ${step4} --rootODir ${rootODir} --release September2015

