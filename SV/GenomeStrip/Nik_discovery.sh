#!/bin/bash

# If you adapt this script for your own use, you will need to set these two variables based on your environment.
# SV_DIR is the installation directory for SVToolkit - it must be an exported environment variable.
# SV_TMPDIR is a directory for writing temp files, which may be large if you have a large data set.
cd /SAN/mottlab/Mitochondria/5_Nik/GenomeStrip
export SV_DIR=`cd /SAN/mottlab/heterosis/bin/svtoolkit && pwd`
SV_TMPDIR=./tmpdir

runDir=chr5_FR07919256
inputFile=/SAN/vyplab/UKIRDC/b37/chr5_FR07919256_sorted_unique.bam
sites=chr5_FR07919256.discovery.vcf
genotypes=chr5_FR07919256.genotypes.vcf

# These executables must be on your path.
which java > /dev/null || exit 1
which Rscript > /dev/null || exit 1
which samtools > /dev/null || exit 1

# For SVAltAlign, you must use the version of bwa compatible with Genome STRiP.
export PATH=${SV_DIR}/bwa:${PATH}
export LD_LIBRARY_PATH=${SV_DIR}/bwa:${LD_LIBRARY_PATH}

mx="-Xmx4g"
classpath="${SV_DIR}/lib/SVToolkit.jar:${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar:${SV_DIR}/lib/gatk/Queue.jar"

mkdir -p ${runDir}/logs || exit 1
mkdir -p ${runDir}/metadata || exit 1

# Display version information.
java -cp ${classpath} ${mx} -jar ${SV_DIR}/lib/SVToolkit.jar

# Run preprocessing.
# For large scale use, you should use -reduceInsertSizeDistributions, but this is too slow for the installation test.
# The method employed by -computeGCProfiles requires a GC mask and is currently only supported for human genomes.
java -cp ${classpath} ${mx} \
    org.broadinstitute.gatk.queue.QCommandLine \
    -S ${SV_DIR}/qscript/SVPreprocess.q \
    -S ${SV_DIR}/qscript/SVQScript.q \
    -gatk ${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar \
    --disableJobReport \
    -cp ${classpath} \
    -configFile ${SV_DIR}/conf/genstrip_parameters.txt \
    -tempDir ${SV_TMPDIR} \
    -R /SAN/vyplab/UKIRDC/reference/human_g1k_v37.fasta \
	-R /SAN/vyplab/UKIRDC/reference/human_g1k_v37.fasta
    -runDirectory ${runDir} \
    -md ${runDir}/metadata \
    -disableGATKTraversal \
    -useMultiStep \
    -reduceInsertSizeDistributions false \
    -computeGCProfiles true \
    -computeReadCounts true \
    -jobLogDir ${runDir}/logs \
    -I ${inputFile} \
    -run \
    || exit 1

# Run discovery.
java -cp ${classpath} ${mx} \
    org.broadinstitute.gatk.queue.QCommandLine \
    -S ${SV_DIR}/qscript/SVDiscovery.q \
    -S ${SV_DIR}/qscript/SVQScript.q \
    -gatk ${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar \
    --disableJobReport \
    -cp ${classpath} \
    -configFile ${SV_DIR}/conf/genstrip_parameters.txt \
    -tempDir ${SV_TMPDIR} \
    -R /SAN/vyplab/UKIRDC/reference/human_g1k_v37.fasta \
    -runDirectory ${runDir} \
    -md ${runDir}/metadata \
    -disableGATKTraversal \
    -jobLogDir ${runDir}/logs \
    -L 1 \
    -minimumSize 100 \
    -maximumSize 1000000 \
    -suppressVCFCommandLines \
    -I ${inputFile} \
    -O ${sites} \
    -run \
    || exit 1

(grep -v ^##fileDate= ${sites} | grep -v ^##source= | grep -v ^##reference= | diff -q - benchmark/${sites}) \
    || { echo "Error: test results do not match benchmark data"; exit 1; }

# Run genotyping on the discovered sites.
java -cp ${classpath} ${mx} \
    org.broadinstitute.gatk.queue.QCommandLine \
    -S ${SV_DIR}/qscript/SVGenotyper.q \
    -S ${SV_DIR}/qscript/SVQScript.q \
    -gatk ${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar \
    --disableJobReport \
    -cp ${classpath} \
    -configFile ${SV_DIR}/conf/genstrip_parameters.txt \
    -tempDir ${SV_TMPDIR} \
    -R /SAN/vyplab/UKIRDC/reference/human_g1k_v37.fasta \
    -runDirectory ${runDir} \
    -md ${runDir}/metadata \
    -disableGATKTraversal \
    -jobLogDir ${runDir}/logs \
    -I ${inputFile} \
    -vcf ${sites} \
    -O ${genotypes} \
    -run \
    || exit 1

(grep -v ^##fileDate= ${genotypes} | grep -v ^##source= | grep -v ^##contig= | grep -v ^##reference= | diff -q - benchmark/${genotypes}) \
    || { echo "Error: test results do not match benchmark data"; exit 1; }

