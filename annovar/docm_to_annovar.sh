# Convert DOCM VCF to ANNOVAR format.

INPUT=$1

HOME_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Create annotation for all mutations.
${HOME_DIR}/vcf_to_annovar.sh $INPUT > hg19_DOCM.txt

# Create subsets of annotations.
selections=("(lung cancer)|(non-small cell.*lung)" "lung adenocarcinoma" "renal")
names=(LUNG_OR_NSCLC LUAD RENAL)
n=${#selections[@]}
for (( i=0; i<${n}; i++ ));
do
    egrep "${selections[$i]}" $INPUT | ${HOME_DIR}/vcf_to_annovar.sh > hg19_DOCM_${names[$i]}.txt
done

