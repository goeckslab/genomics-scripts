#!/bin/bash

#
# Find somatic variants in a VCF using GEMINI.
#
# Dependencies that should be in your PATH:
#   GEMINI (https://github.com/arq5x/gemini, See ./install/install_gemini.sh for installation instructions)
#	vt (https://github.com/atks/vt, see http://genome.sph.umich.edu/wiki/Vt#Installation for installation instructions)
#   bcftools (https://github.com/samtools/bcftools, see http://www.htslib.org/download/ for installation instructions)
#
# Arguments:
#   genome_fasta - fasta file for genome
#   annotator_dir - directory for annotation tool
#   annotations_dir - directory for annotations
#   gemini_pythonpath - PYTHONPATH to GEMINI's python source
#   depth - minimum variant depth [CURRENTLY IGNORED]
#   ab_novel - minimum allelic balance for novel variants. [CURRENTLY IGNORED]
#   counts file - file to write counts data to.
#

# Arguments check.
if [ $# -ne "7" ]
then
  echo "Usage: `basename $0` <genome_fasta> <annotator_dir> <annotations_dir> <gemini_pythonpath> <depth> <ab_novel> <counts_file>"
  exit -1
fi

REFERENCE=$1
ANNOTATOR_DIR=$2
ANNOTATIONS_DIR=$3
GEMINI_PYTHONPATH=$4
DEPTH=$5
AB_NOVEL=$6
COUNTS_FILE=$7

JOBS=2
HOME_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
TEMP="temp.txt"

# Set up directory and file for results.
mkdir find_somatic
pushd find_somatic
touch ${TEMP}

# Find somatic variants for all VCFs.
parallel -j ${JOBS} "${HOME_DIR}/find_somatic_with_gemini.sh {} ${REFERENCE} ${ANNOTATOR_DIR} ${ANNOTATIONS_DIR} ${GEMINI_PYTHONPATH} ${DEPTH} ${AB_NOVEL} ${TEMP}" ::: ../*.vcf

# Create final counts file.
(head -n 1 ${TEMP}; grep -v ^# ${TEMP} | sort -k1,1) > ${COUNTS_FILE}

# Merge all somatic variants into single file.
parallel -j ${JOBS} "bgzip {} && tabix -p vcf {}.gz" ::: *_somatic.vcf
bcftools merge *.vcf.gz > all_somatic.vcf

# Create GEMINI database.
# TODO: annotate with TCGA, DOCM annotations.

# Because of the INFO bug in GEMINI v0.15.1, need to create gemini database (annotation + create)
# rather than just create.
#gemini load -v all_somatic.vcf -t VEP all_somatic.db
${HOME_DIR}/../pipelines/create_gemini_db.sh all_somatic.vcf all_somatic.db VEP ${ANNOTATOR_DIR}

popd