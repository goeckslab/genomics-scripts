#!/bin/sh

# Parameter checking.
if [ $# -ne "3" ]
then
  echo "Usage: `basename $0` <input.vcf> <db_name> <VEP/snpEff>"
  exit -1
fi

# Set up name for annotated VCF.
BASE=$(basename "$1" .vcf)
ANNO_VCF="${BASE}.anno.vcf"
SNPEFF_DIR="/Users/jeremy/tmp/gemini/snpEff"
VEP_DIR="/Users/jeremy/tmp/gemini/vep/ensembl-tools-release-79/scripts/variant_effect_predictor"

# Annotate.
ANNOTATION=$3
if [ $ANNOTATION = "VEP" ]; then
	perl ${VEP_DIR}/variant_effect_predictor.pl -i $1 \
    --cache \
    --offline \
    --assembly GRCh37 \
    --sift b \
    --polyphen b \
    --symbol \
    --numbers \
    --biotype \
    --total_length \
    -o ${ANNO_VCF} \
    --vcf \
    --fields Consequence,Codons,Amino_acids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE
elif [ $ANNOTATION = "snpEff" ]; then
	java -jar ${SNPEFF_DIR}/snpEff.jar -i vcf -o vcf GRCh37.75 $1 > ${ANNO_VCF}
fi

# Load into GEMINI.
gemini load -v ${ANNO_VCF} -t ${ANNOTATION} $2
