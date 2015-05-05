#!/bin/sh

# Parameter checking.
if [ $# -ne "2" ]
then
  echo "Usage: `basename $0` <input.vcf> <db_name>"
  exit -1
fi

# Annotate variants with SnpEff.
BASE=$(basename "$1" .vcf)
java -jar snpeff/snpEff.jar -i vcf -o vcf GRCh37.75 $1 > ${BASE}.snpeff.vcf

# Load into GEMINI.
gemini load -v ${BASE}.snpeff.vcf -t snpEff $2
