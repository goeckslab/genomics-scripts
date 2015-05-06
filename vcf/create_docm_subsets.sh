#
# Create variant subsets of DOCM (http://docm.genome.wustl.edu)
#

# Egrep a VCF, preserving headers.
function egrep_vcf {
	(grep ^# "$2"; grep -v ^# "$2" | egrep "$1")
}

mkdir DOCM
pushd DOCM

# Download DOCM variants.
VARIANTS=DOCM.vcf
wget -O ${VARIANTS}.orig "http://docm.genome.wustl.edu/api/v1/variants.vcf?undefined=&chromosomes=&genes=&diseases=&mutation_types=&amino_acids=&position_start=&position_stop="
(grep ^# ${VARIANTS}.orig; grep -v ^# ${VARIANTS}.orig | sort -k1,1 -k2,2n | sed 's/^/chr/') > ${VARIANTS}

# Create variant subsets.
selections=("(lung cancer)|(non-small cell.*lung)" "lung adenocarcinoma" "lung squamous" "renal")
names=(LUNG_OR_NSCLC LUAD LUSQ RENAL)
n=${#selections[@]}
for (( i=0; i<${n}; i++ )); do
    egrep_vcf "${selections[$i]}" $VARIANTS > "DOCM_${names[$i]}.vcf"
done

# Compress and index.
for f in *.vcf; do
	bgzip $f && tabix -p vcf $f.gz
done

popd

