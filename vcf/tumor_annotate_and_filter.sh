#!/bin/sh
#SBATCH -t 01:30:00

module load parallel

JOBS=16

# 1. Process variants so that there is one variant per line (vcfallelicprimitives, vcfreakmulti), filter to remove variants with reference genotype (these are created when processing variants), 
# and name sample using filename.
parallel -j ${JOBS} "vcfallelicprimitives -k {} | vcfbreakmulti | /groups/cbi/jgoecks/projects/genomics-scripts/vcf/filter_ref_gts.awk | sed 's/unknown/{.}/' > {.}_norm.vcf" ::: *.vcf

# 2. Filter by allelic balance, proximity to homopolymer runs, and depth.
parallel -j ${JOBS} "vcffilter -f 'AB > 0.02' -f 'HP < 5' -f 'DP > 200' {} > {.}_AB_2pct.vcf" ::: *norm.vcf

# 3. Annotate with ANNOVAR and fix ExAC type from String to Float
parallel -j ${JOBS} "/groups/cbi/jgoecks/projects/genomics-scripts/annovar/table_annovar.sh {}" ::: *2pct.vcf
parallel -j ${JOBS} "sed -i.bak 's/exac03\,Number=\.\,Type=String/exac03\,Number=\.\,Type=Float/'" ::: *hg19_multianno.vcf 

# 4. BED annotations with:
#   -Personalized Cancer Therapy mutations (http://pct.mdanderson.org/, https://bitbucket.org/wanding/clinsek#markdown-header-prebuilt-target-site-compilations)
parallel -j ${JOBS} "vcfannotate -b ~/pct.bed -k PCT -d 'NA' {} > {1.}_pctm.vcf" ::: *hg19_multianno.vcf

# 5. Filter to remove common and silent variants.
# TODO: add comment about filter by exonic function or figure out how to get '.' to work for vcffilter
parallel -j ${JOBS} "egrep -v 'ExonicFunc.refGene=(\.|synonymous_SNV)' {} | vcffilter -f '1000g2014oct_all < 0.01' -f 'esp6500si_all < 0.01' -f 'exac03 < 0.01' > {1.}_no_common_or_silent.vcf" ::: *pctm.vcf

# 6. Filter to get TCGA variants.
parallel -j ${JOBS} "vcffilter -f 'TCGA = TRUE' {} > {1.}_tcga.vcf" ::: *no_common_or_silent.vcf

# 7. Filter to get DOCM variants.
parallel -j ${JOBS} "vcffilter -f 'DOCM = TRUE' {} > {1.}_docm.vcf" ::: *no_common_or_silent.vcf

# 8. Filter to get PCT variants.
parallel -j ${JOBS} "vcffilter -f '!( PCT = NA )' {} > {1.}_pct.vcf" ::: *no_common_or_silent.vcf

# 9. Filter to get novel variants.
parallel -j ${JOBS} "vcffilter -f 'AB > 0.1' -f '!( TCGA = TRUE )' -f '!( DOCM = TRUE )' -f 'PCT = NA' -f 'SIFT_pred = D' -f 'Polyphen2_HDIV_pred = D' {} > {1.}_novel.vcf" ::: *no_common_or_silent.vcf

# 10. Merge TCGA (hotspots), DOCM, and novel into single file.
parallel -j ${JOBS} "vcfcombine {1.}_tcga.vcf {1.}_docm.vcf {1.}_novel.vcf > {1.}_somatic.vcf" ::: *no_common_or_silent.vcf

# 11. Count TCGA, DOCM, PCT, novel, and total somatic variants.
parallel "echo -n {} ''; grep -v ^# {} | wc -l" ::: *tcga.vcf | sort -k1,1 > tcga_counts.txt
parallel "echo -n {} ''; grep -v ^# {} | wc -l" ::: *novel.vcf | sort -k1,1 > novel_counts.txt
parallel "echo -n {} ''; grep -v ^# {} | wc -l" ::: *no_common_or_silent_docm.vcf | sort -k1,1 > docm_counts.txt
parallel "echo -n {} ''; grep -v ^# {} | wc -l" ::: *no_common_or_silent_pct.vcf | sort -k1,1 > pct_counts.txt
parallel "echo -n {} ''; grep -v ^# {} | wc -l" ::: *somatic.vcf | sort -k1,1 > somatic_counts.txt
(echo '# Sample TCGA Novel DOCM Total PCT'; paste -d' ' tcga_counts.txt novel_counts.txt docm_counts.txt somatic_counts.txt pct_counts.txt | cut -f1,2,4,6,8,10 -d ' ') | sed -r 's/\_norm([A-Za-z0-9\.\_])+//' > total_counts.txt


