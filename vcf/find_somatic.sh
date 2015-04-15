#!/bin/sh
#SBATCH -t 01:30:00

module load parallel

# Arguments check.
if [ $# -ne "2" ]
then
  echo "Usage: `basename $0` <depth> <ab_novel>"
  exit -1
fi

JOBS=4
DP=$1
AB_NOVEL=$2

# TODO:
# *how to handle variants in COSMIC?
# *need to push all indels through regardless of AC count; currently only indels with AC > 1 are kept for somatic final.

# Setup.
mkdir find_somatic_DP${DP}_AB_NOVEL${AB_NOVEL}
pushd find_somatic_DP${DP}_AB_NOVEL${AB_NOVEL}

# Filter to remove common variants, those with low allele frequency, those near HP runs, and theose with low depth.
parallel -j ${JOBS} "vcffilter -f 'AB > 0.02' -f 'HP < 5' -f 'DP > ${DP}' -f '1000g2014oct_all < 0.01' -f 'esp6500si_all < 0.01' -f 'exac03 < 0.01' {} > {1/.}_AB2pct_no_common.vcf" ::: ../annotated/*hg19_multianno.vcf

# Filter to get variant subsets.
annos=(TCGA DOCM DOCM_LUNG_OR_NSCLC DOCM_LUAD)
for anno in ${annos[@]}; do
    parallel -j ${JOBS} "vcffilter -f '${anno} = TRUE' {} > {.}_${anno}.vcf" ::: *_AB2pct_no_common.vcf
done

# Filter to get novel SNPs. 
parallel -j ${JOBS} "egrep -v 'ExonicFunc.refGene=(\.|synonymous_SNV)' {} | vcffilter -f 'AB > ${AB_NOVEL}' -f '!( TCGA = TRUE )' -f '!( DOCM = TRUE )' -f 'SIFT_pred = D' -f 'Polyphen2_HDIV_pred = D' > {.}_novel_snps.vcf" ::: *_AB2pct_no_common.vcf

# Filter to get novel Indels. 
parallel -j ${JOBS} "vcffilter -f 'AB > ${AB_NOVEL}' -f '!( TCGA = TRUE )' -f '!( DOCM = TRUE )' {} | vcffilter -o -f 'TYPE = ins' -f 'TYPE = del' > {.}_novel_indels.vcf" ::: *_AB2pct_no_common.vcf

# Merge TCGA (hotspots), DOCM lung or NSCLC (hotspots), and novel into single file.
parallel -j ${JOBS} "vcfcombine {.}_TCGA.vcf {.}_DOCM_LUNG_OR_NSCLC.vcf {.}_DOCM_LUAD.vcf {.}_novel_snps.vcf {.}_novel_indels.vcf > {.}_somatic.vcf" ::: *_AB2pct_no_common.vcf

# Count TCGA, DOCM*, novel, and total somatic variants.
sets=(TCGA DOCM DOCM_LUNG_OR_NSCLC DOCM_LUAD novel_snps novel_indels somatic)
for a_set in ${sets[@]}; do
    parallel "echo -n {} ''; grep -v ^# {} | wc -l" ::: *${a_set}.vcf | sort -k1,1 > ${a_set}_counts.txt    
done

(echo '#Sample TCGA Novel_SNP Novel_Indel DOCM DOCM_lung DOCM_luad Total'; paste -d' ' TCGA_counts.txt novel_snps_counts.txt novel_indels_counts.txt DOCM_counts.txt DOCM_LUNG_OR_NSCLC_counts.txt DOCM_LUAD_counts.txt somatic_counts.txt | cut -f1,2,4,6,8,10,12 -d ' ') > total_counts.txt

# NOTE: these next two steps do not work well with Freebayes.

# Combine somatic mutations and fixup.
vcfcombine *somatic.vcf | vcffixup > somatic_combined.vcf

# Filter to get final somatic mutations that (a) appear in more than one sample or (b) appear in hotspots.
vcffilter -o -f 'AC > 1' -f 'DOCM_LUNG_OR_NSCLC = TRUE' -f 'DOCM_LUAD = TRUE' -f 'TCGA = TRUE' somatic_combined.vcf > final_somatic_combined.vcf

# Get final somatic mutations for each sample.
parallel -j ${JOBS} "vcfkeepsamples final_somatic_combined.vcf {= s/\_norm([A-Za-z0-9\.\_])+//  =} | vcffixup | vcffilter -f 'AC > 0' > {.}_final.vcf " ::: *somatic.vcf 

# Get counts for final somatic mutations.
parallel "echo -n {} ''; grep -v ^# {} | wc -l" ::: *somatic_final.vcf | sort -k1,1 > somatic_final_counts.txt

# Get counts for novel SNPs in final somatic mutations.
parallel -j ${JOBS} "(echo -n {} ''; vcffilter -f 'TYPE = snp ' -f '!( TCGA = TRUE )' -f '!( DOCM = TRUE )' -f '!( DOCM_LUNG_OR_NSCLC = TRUE )' -f '!( DOCM_LUAD = TRUE )' {} | grep -v ^# | wc -l)" ::: *somatic_final.vcf | sort -k1,1 > final_novel_snp_counts.txt

# Get counts for novel indels in final somatic mutations.
parallel -j ${JOBS} "(echo -n {} ''; vcffilter -f '!( TYPE = snp )' -f '!( TCGA = TRUE )' -f '!( DOCM = TRUE )' -f '!( DOCM_LUNG_OR_NSCLC = TRUE )' -f '!( DOCM_LUAD = TRUE )' {} | grep -v ^# | wc -l)" ::: *somatic_final.vcf | sort -k1,1 > final_novel_indel_counts.txt

# Get counts for each annotation set for each final somatic mutations.
annos=(TCGA DOCM DOCM_LUNG_OR_NSCLC DOCM_LUAD)
for anno in ${annos[@]}; do
    parallel -j ${JOBS} "echo -n {} ''; vcffilter -f '${anno} = TRUE' {} | grep -v ^# | wc -l" ::: *somatic_final.vcf | sort -k1,1 > ${anno}_final_counts.txt
done

# Merge counts for final somatic mutations.
(echo '#Sample TCGA NovelSNP NovelIndel DOCM DOCM_lung DOCM_luad Total'; paste -d' ' TCGA_final_counts.txt final_novel_snp_counts.txt final_novel_indel_counts.txt DOCM_final_counts.txt DOCM_LUNG_OR_NSCLC_final_counts.txt DOCM_LUAD_final_counts.txt somatic_final_counts.txt | cut -f1,2,4,6,8,10,12 -d ' ') > final_total_counts.txt

# Keep only sample names in count files.
sed -i.bak -r 's/\_norm([A-Za-z0-9\.\_])+//' *.txt
rm *.bak

# Cleanup.
popd