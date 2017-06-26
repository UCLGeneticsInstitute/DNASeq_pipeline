rootODir=${1-$rootODir}
release=${2-$release}
bed=${rootODir}/allChr_snpStats_out.bed
srt=${rootODir}ExomeDepth/allChr_snpStats_sorted.bed
bedtools=/cluster/project8/vyp/cian/support/Bedtools/Bedtools_2.22/bin/mergeBed

oDir=${rootODir}ExomeDepth/
if [ ! -e $oDir ]; then mkdir $oDir; fi

sort -k1,1 -k2,2n $bed > $srt
$bedtools  -i $srt  -d 1000 >  ${oDir}UCLex_merged.bed
