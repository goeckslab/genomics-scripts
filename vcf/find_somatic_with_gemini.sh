#
# Find somatic variants in a VCF using GEMINI.
#
# Dependencies that should be in your PATH:
#   GEMINI
# 	bgzip & tabix
#	vt
#
# Arguments:
#   vcf - VCF file to process
#   genome_fasta - fasta file for genome
#   annotations_dir - directory for annotations
#   gemini_pythonpath - PYTHONPATH to GEMINI's python source
#   depth - minimum variant depth [CURRENTLY IGNORED]
#   ab_novel - minimum allelic balance for novel variants. [CURRENTLY IGNORED]
#   counts file - file to write counts data to.
#

# Arguments check.
if [ $# -ne "7" ]
then
  echo "Usage: `basename $0` <vcf> <genome_fasta> <annotations_dir> <gemini_pythonpath> <depth> <ab_novel> <counts_file>"
  exit -1
fi

# Set up variables.
INPUT_VCF=$1
REFERENCE=$2
ANNOTATIONS_DIR=$3
GEMINI_PYTHONPATH=$4
COUNTS_FILE=$7
HOME_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && cd .. && pwd )"
CREATE_GEMINI_DB="${HOME_DIR}/pipelines/create_gemini_db.sh"
BASE=$(basename "$1" .vcf)
FIND_SOMATIC_SCRIPT="${HOME_DIR}/vcf/gemini_find_somatic.py"

# Remove cruft (low AB, low DP, or high HP), then decompose and normalize.
(grep ^# ${INPUT_VCF}; vt view -f "INFO.AB>0.02&&INFO.DP>50&&INFO.HP<5" ${INPUT_VCF}) | vt decompose -s - | vt normalize -r ${REFERENCE} -o ${BASE}_decnorm.vcf -

# Annotate and create GEMINI db.
${CREATE_GEMINI_DB} ${BASE}_decnorm.vcf ${BASE}.db VEP

# Add annotations.
ANNOS=""
for d in ${ANNOTATIONS_DIR}/*/ ; do
	for f in $d/*.vcf.gz ; do
		ANNO=$(basename "$f" .vcf.gz)
		ANNOS="${ANNOS},${ANNO}"
		gemini annotate \
			-f $f \
			-a boolean \
			-c ${ANNO} \
			${BASE}.db
	done
done

# Create VCF of somatic variants.
PYTHONPATH=${GEMINI_PYTHONPATH} python ${FIND_SOMATIC_SCRIPT} --header ${BASE}.db ${ANNOS:1} ${BASE}_somatic.vcf >> ${COUNTS_FILE}
