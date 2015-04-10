#!/bin/sh

# Script requires BEDtools.

# Arguments check.
if [ $# -ne "2" ]
then
  echo "Usage: `basename $0` <maf_file> <bed_file>"
  exit -1
fi

MAF=$1
BED=$2

# Use awk to create interval (BED-like file from MAF),
# then print MAF header, intersect MAF with BED to get overlaps,
# and remove the initial columns to get back to a MAF.
awk '(NR > 2) {print "chr" $5, $6 -1, $7, $0;}' OFS='\t' $MAF \
| (head -n 2 $MAF; bedtools intersect -a stdin -b $BED | cut -f4-)
