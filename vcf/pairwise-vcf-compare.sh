#!/bin/sh
#SBATCH -t 01:30:00

export PERL5LIB=/groups/cbi/jgoecks/tools/vcftools_0.1.12b/perl
module load parallel

parallel -j 16 '/groups/cbi/jgoecks/tools/vcftools_0.1.12b/bin/vcf-compare {1} {2} | grep ^VN | head -n 1' ::: *.sorted.gz ::: *.sorted.gz > pairwise.txt
