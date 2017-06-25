oFolder=$1
bamList=$2
bamDir=$3

script="/SAN/vyplab/UCLex/scripts/DNASeq_pipeline/ciangene/CNV/ExomeDepth/run_CNVs.sh"
samtools="/share/apps/genomics/samtools-1.3.1/bin/samtools"
Rscript="/cluster/project8/vyp/vincent/Software/R-3.3.0/bin/Rscript"
getGenders="/SAN/vyplab/UCLex/scripts/DNASeq_pipeline/ciangene/CNV/ExomeDepth/makeSampleLists.R"
plot='/SAN/vyplab/UCLex/scripts/DNASeq_pipeline/ciangene/CNV/ExomeDepth/plotExomeDepth_optparse.R'
plotX='/SAN/vyplab/UCLex/scripts/DNASeq_pipeline/ciangene/CNV/ExomeDepth/plotExomeDepth_optparse_X.R'
genes="/SAN/vyplab/UCLex/support/all.genes.bed"

echo "Making outDir"
if [ ! -e $oFolder ]; then mkdir -p $oFolder; fi
if [ ! -e ${oFolder}/ChrDepth/ ]; then mkdir -p ${oFolder}/ChrDepth/; fi
if [ ! -e ${oFolder}/male/ ]; then mkdir -p ${oFolder}/male/; fi
if [ ! -e ${oFolder}/female/ ]; then mkdir -p ${oFolder}/female/; fi

# First I want to work out which samples are male and female so I can process X chromosome appropriately.
echo "Getting read depth"
while read bam; do
	echo $bam
	t=$(basename $bam)
	$samtools idxstats ${bam} > ${oFolder}/ChrDepth/${t}.chromdepths.txt
done < $bamList

echo "Predicting genders from read depth"
$Rscript $getGenders --depthDir ${oFolder}/ChrDepth/ --bamDir $bamDir

# Then do CNV calls
$script ${oFolder}/autosomes/ ${bamList} TRUE
$script ${oFolder}/female/ ${oFolder}/ChrDepth/X/Bams_female FALSE
$script ${oFolder}/male/ ${oFolder}/ChrDepth/X/Bams_male FALSE