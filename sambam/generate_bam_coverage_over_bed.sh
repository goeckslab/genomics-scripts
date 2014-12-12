#!/bin/sh
#
# Computes mean read coverage for a set of BAMs in a directory using regions from an input BED file.
#

# Parameter testing and evaluation.
if [ $# -ne "1" ]
then
  echo "Usage: `basename $0` <input.bed>"
  exit -1
fi

# Get script directory.
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

INPUT_BED=$1

# Process BAMs for coverage in BED regions.
for file in *.bam; do 
	# Generate BAM read coverage over genome.
	$DIR/compute_bam_coverage.sh $file

	# Extract coverage information for BED regions.
	bigWigAverageOverBed -minMax $file.bigwig $INPUT_BED $file.bigwig.bedcoverage

	# Extract mean0, min, max coverage for BED regions.
	awk -v OFS='	' '{if (NR == 1) { split(FILENAME, filename, "."); print filename[1] "_mean0", filename[1] "_min", filename[1] "_max" };print $5, $7, $8}' "$file.bigwig.bedcoverage" > "$file.bigwig.bedcoverage.tab"; 
done

# Generate final file.
(echo 'Region'; cut -f4-6 $INPUT_BED) | paste - *.tab > all_coverage.tabular