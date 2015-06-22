"""
Library for using GEMINI's API for somatic variant calling.
"""

import argparse
import os

from collections import namedtuple
from gemini import GeminiQuery, DefaultRowFormat, VCFRowFormat

COMMON_DATABASES = ["1kg", "exac", "esp"]
BASE_VARIANT_QUERY = "select chrom, start, end, ref, alt from variants"

class Variant(object):
    def __init__(self, chrom=None, start=None, end=None, ref=None, alt=None, transcript=None,
                 codon_change=None):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.ref = ref
        self.alt = alt
        self.transcript = transcript
        self.codon_change = codon_change

    def is_adjacent(self, variant):
        """
        Returns true if this variant is next to another variant.
        """
        return self.chrom == variant.chrom and (self.start == variant.end or self.end == variant.start)

    def combine_adjacent(self, variant):
        """
        Combine adjacent variants into a single variant.
        """
        if not self._is_adjecent(variant):
            return None

        if self.start < variant.start:
            first = self
            second = variant
        else:
            first = variant
            second = self

        return Variant(chrom=self.chrom, start=first.start, end=second.end, ref=(first.ref + second.ref),
                       alt=(first.alt + second.alt))

def get_gt_filter(sample, gt_type):
    """
    Returns clause for filtering variants based on sample and GT type.
    """
    return "gt_types.%s == %s" % (sample, gt_type)

def get_gt_count_filter(gt_type, count):
    """
    Returns clause for filtering genotypes based on type and count.
    """
    return "(gt_types).(*).(==%s).(count > %i)" % (gt_type, count)

def get_in_and_aff_clause(db, aaf=0.01):
    """
    Returns clause for checking if variant is below a certain alternate allele frequency in a database
    either because the variant is not in the database or its AAF is below the threshold.
    """
    return "NOT in_%s OR aaf_%s_all < %0.3f" % (db, db, aaf)

def get_no_common_vars_clause(aaf=0.01):
    """
    Returns clause for removing common variants from a query.
    """
    return ' AND '.join( ['(%s)' % get_in_and_aff_clause(db, aaf) for db in COMMON_DATABASES ] )

def get_annotation_clause(anno_col_name, has_anno=True):
    val = 1
    if not has_anno:
        val = 0
    return "%s = %i" % (anno_col_name, val)

def get_annotation_and_no_common_clause(anno_col_name, has_anno=True, aaf=0.01):
    """
    Returns clause for selecting variants with or without an annotation 
    but not appearing in common databases at a given alternate allele frequency.
    """
    val = 1
    if not has_anno:
        val = 0
    return '%s AND %s' % ( get_annotation_clause(anno_col_name, val), get_no_common_vars_clause() )

def get_novel_query(annotations, var_type='snp', allele_freq=0.1):
    """
    Return query to get novel variants.
    """

    add_clause = ""
    if var_type is 'snp':
        #TODO: be more precise about which impacts in HIGH (e.g., stop_gain, stop_loss) and MED (e.g. non-synonymous) 
        # should be selected
        # Add clause for novel SNPs to require high impact or med impact + SIFT/PolyPhen to be true.
        add_clause = "AND (impact_severity = 'HIGH' or (impact_severity = 'MED' AND " \
                     "sift_pred = 'deleterious' AND polyphen_pred = 'probably_damaging'))"
    
    # Set up query.
    no_anno = [get_annotation_clause(anno, has_anno=False) for anno in annotations]
    return "%s WHERE (type = '%s') AND (allele_bal >= %f) AND %s AND %s %s" \
           % ( BASE_VARIANT_QUERY, var_type, allele_freq, " AND ".join(no_anno), get_no_common_vars_clause(), add_clause)

def get_hotspot_variants(annotation, allele_bal=0.02):
    """
    Returns query to select variants in a hotspot, i.e. those with an annotation and minimum allele_bal
    """
    return "%s WHERE (allele_bal >= %f) AND %s" % \
          ( BASE_VARIANT_QUERY, allele_bal, get_annotation_and_no_common_clause(annotation) )

def get_somatic_variants(gemini_db, annotations, out_format=DefaultRowFormat(None)):
    """
    Returns counts for all hotspot and novel variants and a list of variants.
    """

    variants = []
    var_counts = []

    def add_variants(results):
        count = 0
        for result in results:
            result = str(result)
            if result not in variants:
                variants.append(result)
            else:
                pass
                #print "******************* Result already in variants:"
                #print result.split('\t')[0:3]
            count += 1
        return count

    # Get variants in hotspots.
    for anno in annotations:
        results = get_query_results(gemini_db, get_hotspot_variants(anno), out_format=out_format)
        count = add_variants(results)
        var_counts.append(count)

    # Get novel counts.
    for var_type in ['snp', 'indel']:
        results = get_query_results(gemini_db, get_novel_query(annotations, var_type), out_format=out_format)
        count = add_variants(results)
        var_counts.append(count)
    return var_counts, variants

def get_amplicons(gemini_db):
    """
    Returns a tuple of amplicon, variants for all samples and amplicons that share more than one variant.
    """

    query = "SELECT gene, amplicon, gts.%s FROM variants"
    gt_filter = "(gt_types).(*).(==HET).(count > 1) and gt_types.%s == HET or gt_types.%s == HOM_ALT"

    for sample in get_samples(gemini_db):
        variants = get_query_results(gemini_db, query % sample, gt_filter % (sample, sample))
        amplicons_count = {}
        for variant in variants:
            #print variant
            amplicons = variant['amplicon'].split(',')
            for amplicon in amplicons:
                if amplicon not in amplicons_count:
                    amplicons_count[amplicon] = 0
                amplicons_count[amplicon] += 1

        amplicons_with_multiple_variants = False
        for amplicon, count in amplicons_count.items():
            if count > 1:
                print sample, amplicon, count
                amplicons_with_multiple_variants = True

        if not amplicons_with_multiple_variants:
            print sample, 'None', 0

def get_samples(gemini_db):
    """
    Returns list of samples in a GEMINI database.
    """
    return [str(sample['name']) for sample in get_query_results(gemini_db, "select name from samples")]

def has_sample(gemini_db, sample):
    return int( str( get_query_results(gemini_db, "select count(*) from samples where name='%s'" % sample).next() ) ) != 0

def compare_replicates(gemini_db):
    """
    Compare all replicates for shared variants and variants likely caused by deamination.
    """

    query = "SELECT chrom, start, ref, alt, gene, %s, %s FROM variants" 
    gt_filter = "gt_types.%s == HET or gt_types.%s == HET or gt_types.%s == HOM_ALT or gt_types.%s == HOM_ALT"

    # Loop through samples, querying replicates.
    samples = get_samples(gemini_db)
    for sample in samples:
        # Only compare repeats.
        if '_Repeats_' not in sample:
            continue

        # Query for variants where one or both samples are alternate.
        sample_repeat = sample
        sample_original = sample_repeat.replace('_Repeats_', '_FirstBatch_')
        sample_repeat_gts = "gts." + sample_repeat
        sample_original_gts = "gts." + sample_original

        # If original sample is not present, skip.
        if sample_original not in samples:
            #print 'No original:', sample_original
            continue

        #print sample_repeat, sample_original
        #print query % (sample_original_gts, sample_repeat_gts), gt_filter % (sample_original, sample_repeat, sample_original, sample_repeat)

        variants = get_query_results(gemini_db, query % (sample_original_gts, sample_repeat_gts),
                                     gt_filter % (sample_original, sample_repeat, sample_original, sample_repeat))

        shared_count = 0
        unique_count = 0
        deamination_count = 0
        for variant in variants:
            #print variant
            if variant[sample_original_gts] == variant[sample_repeat_gts]:
                shared_count += 1
            elif ( variant['ref'] == 'C' and (variant[sample_original_gts] == 'C/T' or variant[sample_repeat_gts] == 'C/T') ) or \
                 ( variant['ref'] == 'G' and (variant[sample_original_gts] == 'G/A' or variant[sample_repeat_gts] == 'G/A') ):
                deamination_count += 1
            else:
                # TODO: go back to original sample databases and see if variant exists at some AF (e.g. 5%); if so,
                # keep variant. Example: 'chr7  116435999' for NATCH_FirstBatch_MRO10434-0022-M01-36124R3 where
                # FirstBatch AF = 0.06 and Repeats AF = 0.12
                #print variant
                unique_count += 1

        print sample_original, shared_count, deamination_count, unique_count

def query_sample_het(gemini_db, sample, cols="chrom, start, end, ref, alt, gene, cosmic_ids", min_het_count=0, addl_gt_filter=None):
    """
    Query database, returning columns + sample genotype for variants that are (a) HET for sample; 
    (b) have a minimum number of hets total; and (c) meet additional genotype filters.
    """
    query = "select %s, gts.%s from variants" % (cols, sample)
    gt_filter = "gt_types.%s == HET and (gt_types).(*).(==HET).(count >= %i)" % (sample, min_het_count)
    if addl_gt_filter:
        gt_filter += ' and %s' % addl_gt_filter
    return get_query_results(gemini_db, query, gt_filter)

def get_query_results(gemini_db, query, gt_filter="", out_format=DefaultRowFormat(None)):
    """
    Returns results of query.
    """

    gemini = GeminiQuery(gemini_db, out_format=out_format)
    gemini.run(query, gt_filter=gt_filter)
    return gemini

if __name__ == "__main__":
    # Argument setup and parsing.
    parser = argparse.ArgumentParser()
    parser.add_argument("operation", help="Operation to perform")
    parser.add_argument("gemini_db", help="Gemini database to use")
    parser.add_argument("--sample", help="Sample to query for")
    parser.add_argument("--cols", help="Columns to query for")
    parser.add_argument("--gt_count", help="Minimum HET count")
    parser.add_argument("--annotations", help="Annotations to query for")
    parser.add_argument("--output_vcf", help="Write variants to this file")
    parser.add_argument("--header", action="store_true", help="Print header?")
    args = parser.parse_args()
    operation = args.operation
    gemini_db = args.gemini_db
    
    # Do operation.
    if operation == "find_somatic":
        annotations = args.annotations.split(',')
        db_name, ext = os.path.splitext( os.path.split(gemini_db)[1] )
        output = open(args.output_vcf, 'w')

        # Set up VCF formatter.
        simple_struct = namedtuple('Simple', 'db')
        vcf_format = VCFRowFormat(simple_struct(db=gemini_db))
        output.write(vcf_format.header(None) + '\n')

        # Get counts, variants.
        counts, variants = get_somatic_variants(gemini_db, annotations, out_format=vcf_format)

        # Sort variants by chrom and start position.
        def get_chrom_and_start(v):
            fields = v.split('\t')
            return (int(fields[0][3:]), int(fields[1]))
        variants.sort(key=get_chrom_and_start)
        
        # Print header?
        if args.header:
            print '#Sample %s Novel_SNPS Novel_Indels' % ' '.join(annotations)

        # Print output of variant counts.
        print '%s %s' % ( db_name, ' '.join(str(c) for c in counts) )

        # Print somatic variants VCF.
        for variant in variants:
            output.write(variant + '\n')

        output.close()
    elif operation == "compare_replicates":
        print "TODO: connect to function"
    elif operation == "print_samples":
        for sample in get_samples(gemini_db):
            print sample
    elif operation == "query_sample":
        sample = args.sample
        cols = args.cols or "chrom, start, end, ref, alt, gene, cosmic_ids"
        gt_count = args.gt_count or 0
        gt_count = int(gt_count)
        for row in query_sample_het(gemini_db, sample, cols, gt_count):
            print row


