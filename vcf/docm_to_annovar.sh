# Convert DOCM VCF to ANNOVAR format.

INPUT=$1

HOME_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Create annotation for all mutations.
${HOME_DIR}/vcf_to_annovar.sh $INPUT > hg19_DOCM.txt

# Create annotation for Lung or NSCLC mutations.
egrep '(lung cancer)|(non-small cell.*lung)' variants.vcf | ${HOME_DIR}/vcf_to_annovar.sh > hg19_DOCM_LUNG_OR_NSCLC.txt

# Create annotation for lung adenocarcinoma.
egrep 'lung adenocarcinoma' variants.vcf | ${HOME_DIR}/vcf_to_annovar.sh > hg19_DOCM_LUAD.txt
