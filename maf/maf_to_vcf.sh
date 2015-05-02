# Arguments check.
if [ $# -ne "2" ]
then
  echo "Usage: `basename $0` <input_maf> <output_vcf>"
  exit -1
fi

# Add header.
echo '##fileformat=VCFv4.1' > $2
echo '#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO' >> $2

# Convert MAF to VCF by printing VCF columns and sorting and removing duplicates.
tail -n +3 $1 | awk -v FS='	' -v OFS='	' "{print \"chr\"\$5, \$6, \".\", \$11, \$13, \".\", \".\", \"TCGA_LUAD=TRUE\"}" | sort -k1,1 -k2,2n | uniq >> $2
