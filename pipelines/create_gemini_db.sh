#!/bin/sh

#
# Create a GEMINI database by (a) annotating variants with VEP or snpEff and (b) loading annotated variants
# into GEMINI.
#
# Dependencies:
#   VEP and/or snpEff
#   gemini
#

# Parameter checking.
if [ $# -ne "4" ]
then
  echo "Usage: `basename $0` <input.vcf> <db_name> <VEP/snpEff> <VEP/snpEff directory>"
  exit -1
fi

# Set up name for annotated VCF.
BASE=$(basename "$1" .vcf)
ANNO_VCF="${BASE}.anno.vcf"
GEMINI_DB=$2
ANNOTATION=$3
ANNO_DIR=$4

# Annotate.
if [ $ANNOTATION = "VEP" ]; then
    perl ${ANNO_DIR}/variant_effect_predictor.pl -i $1 \
    --cache \
    --refseq \
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
	java -jar ${ANNO_DIR}/snpEff.jar -i vcf -o vcf GRCh37.75 $1 > ${ANNO_VCF}
fi

# Load into GEMINI.
gemini load -v ${ANNO_VCF} -t ${ANNOTATION} ${GEMINI_DB}
