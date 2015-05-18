#!/bin/bash

#
# Find somatic variants in a VCF using GEMINI.
#
# Dependencies that should be in your PATH:
#   GEMINI
#	vt
#
# Arguments:
#   genome_fasta - fasta file for genome
#   annotations_dir - directory for annotations
#   gemini_pythonpath - PYTHONPATH to GEMINI's python source
#   depth - minimum variant depth [CURRENTLY IGNORED]
#   ab_novel - minimum allelic balance for novel variants. [CURRENTLY IGNORED]
#   counts file - file to write counts data to.
#

# Arguments check.
if [ $# -ne "6" ]
then
  echo "Usage: `basename $0` <genome_fasta> <annotations_dir> <gemini_pythonpath> <depth> <ab_novel> <counts_file>"
  exit -1
fi

REFERENCE=$1
ANNOTATIONS_DIR=$2
GEMINI_PYTHONPATH=$3
DEPTH=$4
AB_NOVEL=$5
COUNTS_FILE=$6

JOBS=2
HOME_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
TEMP="temp.txt"

# Set up directory and file for results.
mkdir find_somatic
pushd find_somatic
touch ${TEMP}

# Find somatic variants for all VCFs.
parallel -j ${JOBS} "${HOME_DIR}/find_somatic_with_gemini.sh {} ${REFERENCE} ${ANNOTATIONS_DIR} ${GEMINI_PYTHONPATH} ${DEPTH} ${AB_NOVEL} ${TEMP}" ::: ../*.vcf

# Create final counts file.
(head -n 1 ${TEMP}; grep -v ^# ${TEMP} | sort -k1,1) > ${COUNTS_FILE}

# Merge all somatic variants into single file.
parallel -j ${JOBS} "bgzip {} && tabix -p vcf {}.gz" ::: *_somatic.vcf
bcftools merge *.vcf.gz > all_somatic.vcf

# Create GEMINI database.
# TODO: annotate with TCGA, DOCM annotations.
gemini load -v all_somatic.vcf -t VEP all_somatic.db

popd