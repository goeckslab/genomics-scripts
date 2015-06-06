#!/bin/bash
#
# Compute read coverage for a BAM.
#
# Usage:
#   compute_bam_coverage.sh <input.bam>
#
# Dependencies: bedtools

# Parameter testing and evaluation.
if [ $# -ne "1" ]
then
  echo "Usage: `basename $0` <input.bam>"
  exit -1
fi

# Get directory of script, which is used for finding bam_to_chrom_len
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Create chrom len file.
$DIR/bam_to_chrom_len.sh $1 > genome.len

# Compute coverage. 
bedtools genomecov -bg -split -ibam $1 -g genome.len | wigToBigWig stdin genome.len $1.bigwig

# Clean up.
rm genome.len