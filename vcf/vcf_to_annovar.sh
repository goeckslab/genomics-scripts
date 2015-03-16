# Convert a VCF into a file for use with ANNOVAR.

INPUT=${1:-/dev/stdin}
grep -v ^# ${INPUT} | awk -v OFS='\t' '{ print "chr" $1, $2, $2 + length($5) - 1, $4, $5, "TRUE"}'
