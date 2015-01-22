#!/bin/bash

# Setup
module load parallel
export PERL5LIB=/groups/cbi/jgoecks/tools/vcftools_0.1.12b/perl
# Change this to 'cat' to avoid preprocessing.
PREPROCESS=/groups/cbi/jgoecks/tools/bin/vcfallelicprimitives

# 1. Copy sample and repeat to directory.

# 2. Preprocess VCF, e.g. to create primitives.
parallel "${PREPROCESS} {} > {.}_prim.vcf" ::: *.vcf

# 3. Filter vcfs by depth and quality:
parallel "vcffilter -f 'DP > {2}' -f 'QUAL > {3}' {1} > {1.}_DP{2}_QUAL{3}.vcf" ::: *_prim.vcf ::: 0 20 50 100 ::: 0 1 5 10 20

# 4. Sort, compress, and index all VCF files in a directory.
parallel '(grep ^# {}; grep -v ^# {} | sort -k1,1 -k2,2n) | bgzip > {}.sorted.gz; tabix -p vcf {}.sorted.gz' ::: *_prim*.vcf

# 5. Generate comparisons between first batch and repeats.
parallel -j 4 '/groups/cbi/jgoecks/tools/vcftools_0.1.12b/bin/vcf-compare {1} {2} | grep ^VN' ::: *FirstBatch*_prim*DP*.sorted.gz ::: *Repeats*_prim*DP*.sorted.gz > pairwise.txt

# 5. Sort and view.
#sort -k2,2n pairwise.txt | less
#sort -k2,2n pairwise.txt | egrep 'DP20_.*DP20_'