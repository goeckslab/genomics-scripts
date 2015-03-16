#!/bin/sh
#SBATCH -t 01:30:00

module load parallel

JOBS=16
DP=200

#TODO: make AB, DP parameter configurable, try 1% for hotspots.

# Process variants so that there is one variant per line (vcfallelicprimitives, vcfreakmulti), filter to remove variants with reference genotype (these are created when processing variants), 
# and name sample using filename.
parallel -j ${JOBS} "vcfallelicprimitives -k {} | vcfbreakmulti | sed 's/unknown/{.}/' > {.}_norm.vcf" ::: *.vcf

# Filter by allelic balance, proximity to homopolymer runs, and depth.
parallel -j ${JOBS} "vcffilter -f 'AB > 0.02' -f 'HP < 5' -f 'DP > ${DP}' {} > {.}_DP_${DP}_AB_2pct.vcf" ::: *norm.vcf

# Annotate with ANNOVAR and fix ExAC type from String to Float
parallel -j ${JOBS} "/groups/cbi/jgoecks/projects/genomics-scripts/annovar/table_annovar.sh {}" ::: *2pct.vcf
parallel -j ${JOBS} "sed -i.bak 's/exac03\,Number=\.\,Type=String/exac03\,Number=\.\,Type=Float/'" ::: *hg19_multianno.vcf 

# Filter to remove common and silent variants.
# TODO: add comment about filter by exonic function or figure out how to get '.' to work for vcffilter
parallel -j ${JOBS} "egrep -v 'ExonicFunc.refGene=(\.|synonymous_SNV)' {} | vcffilter -f '1000g2014oct_all < 0.01' -f 'esp6500si_all < 0.01' -f 'exac03 < 0.01' > {1.}_no_common_or_silent.vcf" ::: *hg19_multianno.vcf

# Filter to get TCGA variants.
parallel -j ${JOBS} "vcffilter -f 'TCGA = TRUE' {} > {1.}_tcga.vcf" ::: *no_common_or_silent.vcf

# Filter to get DOCM variants.
parallel -j ${JOBS} "vcffilter -f 'DOCM = TRUE' {} > {1.}_docm.vcf" ::: *no_common_or_silent.vcf

# Filter to get DOCM_LUNG_OR_NSCLC variants.
parallel -j ${JOBS} "vcffilter -f 'DOCM_LUNG_OR_NSCLC = TRUE' {} > {1.}_docm_lung_or_nsclc.vcf" ::: *no_common_or_silent.vcf

# Filter to get DOCM_LUNG_LUAD variants.
parallel -j ${JOBS} "vcffilter -f 'DOCM_LUAD = TRUE' {} > {1.}_docm_luad.vcf" ::: *no_common_or_silent.vcf

# Filter to get novel variants.
parallel -j ${JOBS} "vcffilter -f 'AB > 0.1' -f '!( TCGA = TRUE )' -f '!( DOCM = TRUE )' -f 'SIFT_pred = D' -f 'Polyphen2_HDIV_pred = D' {} > {1.}_novel.vcf" ::: *no_common_or_silent.vcf

# Merge TCGA (hotspots), DOCM lung or NSCLC (hotspots), and novel into single file.
parallel -j ${JOBS} "vcfcombine {1.}_tcga.vcf {1.}_docm_lung_or_nsclc.vcf {1.}_novel.vcf > {1.}_somatic.vcf" ::: *no_common_or_silent.vcf

# Count TCGA, DOCM*, novel, and total somatic variants.
parallel "echo -n {} ''; grep -v ^# {} | wc -l" ::: *tcga.vcf | sort -k1,1 > tcga_counts.txt
parallel "echo -n {} ''; grep -v ^# {} | wc -l" ::: *novel.vcf | sort -k1,1 > novel_counts.txt
parallel "echo -n {} ''; grep -v ^# {} | wc -l" ::: *docm.vcf | sort -k1,1 > docm_counts.txt
parallel "echo -n {} ''; grep -v ^# {} | wc -l" ::: *docm_lung_or_nsclc.vcf | sort -k1,1 > docm_lung_or_nsclc_counts.txt
parallel "echo -n {} ''; grep -v ^# {} | wc -l" ::: *docm_luad.vcf | sort -k1,1 > docm_luad_counts.txt
parallel "echo -n {} ''; grep -v ^# {} | wc -l" ::: *somatic.vcf | sort -k1,1 > somatic_counts.txt
(echo '#Sample TCGA Novel DOCM DOCM_lung DOCM_luad Total'; paste -d' ' tcga_counts.txt novel_counts.txt docm_counts.txt docm_lung_or_nsclc_counts.txt docm_luad_counts.txt somatic_counts.txt | cut -f1,2,4,6,8,10,12 -d ' ') | sed -r 's/\_norm([A-Za-z0-9\.\_])+//' > total_counts.txt


