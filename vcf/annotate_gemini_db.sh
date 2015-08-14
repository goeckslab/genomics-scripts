# Add annotations to a GEMINI database. Annotations are in a directory and must
# be in VCF format, bgzip'ed, and tabix'ed.

# Arguments check.
if [ $# -lt "1" ]
then
  echo "Usage: `basename $0` <gemini_db> [annotations_dir] "
  exit -1
fi

# Get script directory.
HOME_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Get shared settings.
source ${HOME_DIR}/settings.sh

# Set inputs and parameters.
GEMINI_DB=$1
ANNOTATIONS_DIR=${2:-${SHARED_ANNOTATIONS_DIR}}

# Annotate.
for d in ${ANNOTATIONS_DIR}/*/ ; do
	for f in $d/*.vcf.gz ; do
		ANNO=$(basename "$f" .vcf.gz)
		gemini annotate \
			-f $f \
			-a boolean \
			-c ${ANNO} \
			${GEMINI_DB}
	done
done
