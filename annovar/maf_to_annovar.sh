#!/bin/bash

# Usage:
#   maf_to_annovar.sh <input.maf> > hg19_TCGA_RCC.txt
#   maf_to_annovar.sh < <input.maf> >  hg19_TCGA_RCC.txt


# Prepare a MAF file (e.g., from TCGA) for use with ANNOVAR.
# tail +3 to skip two line header, then cut chrom, start, end, ref allele, and tumor alt allele. Finally sort and keep unique.
tail -n +3 ${1:-/dev/stdin} | cut -f5-7,11,13 | awk -v OFS='\t' '($4 != $5) { print "chr"$0, "TRUE" }' | sort -k1,1n -k2,2 | uniq
