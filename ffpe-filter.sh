#!/bin/bash
# Filter FFPE artifacts from MAF file of called somatic SNVs

set -euo pipefail


# DEPENDENCIES ###############################################################

# NOTE: Edit the below paths according to your environment.

# Matlab (2013a) Common Runtime installation directory
#mcr_root=/opt/mcr/v81
mcr_root=/broad/software/nonfree/Linux/redhat_5_x86_64/pkgs/matlab_2013a

# Java 1.7
java=java

# Python 2.7
# Ensure that `python` is is available on $PATH


# INPUTS #####################################################################

if [[ $# < 6 ]]; then
  echo "usage: ${##*/} <id> <bam> <maf> <ref> <dbsnp> <outdir>"
  exit 1
fi

id=$1            # sample id
bam=$2           # path to the BAM file
maf=$3           # path to the MAF file containing called variants
ref=$4           # path to the reference fasta
dbsnp=$5         # path to the dbSNP reference
outdir=$6        # output directory


# PREAMBLE ###################################################################

mkdir -p $outdir

# Directory in which this script resides
root="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"


# ANALYSIS ###################################################################

# Task 1
# Calculate the overall FFPE quality score across all reads

if [[ ! -f "$outdir/${id}.ffpe_metrics.txt" ]]; then
  $root/script/1_annotateFFPEbiasQ/annotateFFPEbias.sh \
     -J $java \
     -j $root/jar/CollectOxoGMetrics.FFPE.jar \
     -i ${id} \
     -b ${bam} \
     -r ${ref} \
     -d ${dnsnp} \
     -f G \
     -t A \
     -c CG. \
     -x ffpe_metrics \
     -s 40000000 \
     -o $outdir \

fi


# Task 3
# Append the FFPE quality score to the MAF file

value=$(cat "$outdir/${id}.ffpe_metrics.txt")

$root/script/2_AppendAnnotation2MAF/AppendAnnotation2MAF.sh \
   -i ${id} \
   -m ${maf} \
   -f ffpe_Q \
   -v $value \
   -o $outdir \


# Task 4
# Count reads in the F1R2 and F2R1 configurations
# supporting the ref and alt alleles

inputFile=$outdir/${id}.ffpe_Q.maf.annotated
outputFilename=${id}.ffpe_Q.maf.filtered
poxoG=0.96
artifactThresholdRate=0.01
logThresholdRate=-1
biasQP1=30
biasQP2=1.5
reference_allele=G
artifact_allele=A
bias_field=i_ffpe

$root/script/3_OrientationBias_filter/build/run_orientationBiasFilter.sh $mcr_root \
  ${inputFile} ${outputFilename} $outdir \
  0 0 \
  ${poxoG} \
  ${artifactThresholdRate} \
  ${logThresholdRate} \
  ${biasQP1} \
  ${biasQP2} \
  ${reference_allele} \
  ${artifact_allele} \
  ${bias_field} \

# Final output is:
# ${id}.ffpe_Q.maf.filtered

