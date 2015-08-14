#
# Compress and index all targets.
# 
# Usage: bgzip_and_tabix.sh [directory] [targets]
#
# Default paremeters are current directory (".") and all VCF files ("*.vcf")
#


DIR="${1:-.}"
TARGETS="${2:-*.vcf}"
JOBS=4

pushd ${DIR}
# For sorting:
#parallel -j ${JOBS} '(grep ^# {}; grep -v ^# {} | sort -k1,1 -k2,2n) | bgzip > {}.sorted.gz; tabix -p vcf {}.sorted.gz' ::: ${TARGETS}
parallel -j ${JOBS} 'bgzip {} && tabix -p vcf {}.gz' ::: ${TARGETS}
popd