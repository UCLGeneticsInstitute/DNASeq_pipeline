oFolder=$1
BAMlist=$2
autosome=$3

for folder in ${oFolder}/temp ${oFolder}/cluster ${oFolder}/cluster/R ${oFolder}/cluster/out ${oFolder}/cluster/error ${oFolder}/cluster/submission
do
    if [ ! -e $folder ]; then mkdir -p $folder; fi
done

includeChr=FALSE
iFolder=/SAN/vyplab/UCLex_raw/
everted=FALSE
bedFile=/cluster/project8/vyp/exomeDepth/data/exome_depth_hg19.bed
nCNVs=20

code=`basename $PWD`
script=${oFolder}/cluster/submission/CNVcallsv2_${code}.sh  ##default which can change

until [ -z "$1" ]; do
	# use a case statement to test vars. we always test $1 and shift at the end of the for block.
    case $1 in
	--nCNVs )
	    shift
	    nCNVs=$1;;
	--oFolder )
	    shift
	    oFolder=$1;;
	--code )
	    shift
	    code=$1;;
	--iFolder )
	    shift
	    iFolder=$1;;
	--bedFile )
	    shift
	    bedFile=$1;;
	--script )
	    shift
	    script=$1;;
	--everted )
	    shift
	    everted=$1;;
	--BAMlist)
	    shift
	    BAMlist=$1;;
	-* )
	    echo "Unrecognized option: $1"
	    exit 1;;
    esac
    shift
    if [ "$#" = "0" ]; then break; fi
done

newCode=`basename $script .sh`

if [[ "$BAMlist"  != "none" ]]; then
    BAMoption="--BAMlist=${BAMlist}"
else
    BAMoption=""
fi

if [[ "$bedFile" != "NULL" ]]; then
    BEDoption="--bedFile=${bedFile}"
else
    BEDoption=""
fi


countScript=/SAN/vyplab/UCLex/scripts/DNASeq_pipeline/ciangene/CNV/ExomeDepth/CNV_pipeline/get_depth_v2.R
evertedScript=/SAN/vyplab/UCLex/scripts/DNASeq_pipeline/ciangene/CNV/ExomeDepth/CNV_pipeline/get_everted.R
callScript=/SAN/vyplab/UCLex/scripts/DNASeq_pipeline/ciangene/CNV/ExomeDepth/CNV_pipeline/calling_v2.R
XcallScript=/SAN/vyplab/UCLex/scripts/DNASeq_pipeline/ciangene/CNV/ExomeDepth/CNV_pipeline/call_X.R
Rbin=/cluster/project8/vyp/vincent/Software/R-3.3.0/bin/R

echo "
#!/bin/bash
#$ -S /bin/bash
#$ -o ${oFolder}/cluster/out
#$ -e ${oFolder}/cluster/error
#$ -cwd
#$ -pe smp 1
#$ -l tmem=5.9G,h_vmem=5.9G
#$ -l h_rt=60:0:0
#$ -V
#$ -R y


$Rbin CMD BATCH --no-save --no-restore  --includeChr=${includeChr} ${BEDoption} ${BAMoption} --iFolder=${iFolder} --oFolder=${oFolder} $countScript ${oFolder}/cluster/R/get_count_${newCode}.Rout
" > $script

if [[ "$everted" == "TRUE" ]]; then
echo "
$Rbin CMD BATCH --no-save --no-restore  --oFolder=${oFolder}_everted ${evertedScript} ${oFolder}/cluster/R/get_everted_${newCode}.Rout

" >> $script
fi

if [[ "$autosome" == "TRUE" ]]; then
echo "
$Rbin CMD BATCH --no-save --no-restore --nCNVs=${nCNVs} --oFolder=${oFolder} ${callScript} ${oFolder}/cluster/R/get_calls_${newCode}.Rout

" >> $script
fi

if [[ "$autosome" == "FALSE" ]]; then
echo "
$Rbin CMD BATCH --no-save --no-restore --nCNVs=${nCNVs} --oFolder=${oFolder} ${XcallScript} ${oFolder}/cluster/R/get_calls_${newCode}_X.Rout
" >> $script
fi





echo $script >> runMe