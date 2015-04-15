#!/bin/sh
#SBATCH -t 01:30:00

module load parallel

JOBS=16
HOME_DIR=/groups/cbi/jgoecks/projects/genomics-scripts

mkdir annotated
pushd annotated

# Process variants so that there is one variant per line (vcfallelicprimitives, vcfreakmulti),
# and name sample using filename.
parallel -j ${JOBS} "vcfallelicprimitives -k {} | vcfbreakmulti | sed \"s/unknown/{/.}/\" > {/.}_norm.vcf" ::: ../*.vcf

# Annotate with ANNOVAR and fix ExAC type from String to Float
parallel -j ${JOBS} "${HOME_DIR}/annovar/table_annovar.sh {}" ::: *norm.vcf
parallel -j ${JOBS} "sed -i.bak 's/exac03\,Number=\.\,Type=String/exac03\,Number=\.\,Type=Float/'" ::: *hg19_multianno.vcf 

popd