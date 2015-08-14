#
# Statistic comparison between VCF replicates and non-replicates using a directory of VCF files. Replicates are identified using the following convention:
# '*FirstBatch*' is replicate one and "*Repeats* is replicate two.
#

TARGET_DIRECTORY=${1:-.}
HOME_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && cd ../.. && pwd )"

# Compress and index VCFs.
${HOME_DIR}/vcf/bgzip_and_tabix.sh 

# Merge VCFs.
bcftools merge *.vcf.gz > all.vcf

# Create GEMINI database.
${HOME_DIR}/pipelines/create_gemini_db.sh all.vcf all.db VEP /groups/cbi/jgoecks/tools/ensembl-tools-release-80/scripts/variant_effect_predictor

# Query GEMINI database for genotypes.
gemini query -q "select variant_id, chrom, start, (gts).(*) from variants" --header all.db | \

# Run comparison.
python ${HOME_DIR}/vcf/replicates/replicates_vs_nonreplicates.py
