from __future__ import print_function
import argparse
import os

from collections import namedtuple
from gemini import GeminiQuery, DefaultRowFormat, VCFRowFormat

"""
Utility functions for working with GEMINI's API.
"""

COMMON_DATABASES = ["1kg", "exac", "esp"]
BASE_VARIANT_QUERY = "select chrom, start, end, ref, alt from variants"

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

    # Add clause for novel SNPs to require SIFT/PolyPhen to be true.
    add_clause = None
    if var_type is 'snp':
        add_clause = ""
        add_clause = "AND (sift_pred = 'deleterious' AND polyphen_pred = 'probably_damaging')"
    else:
        add_clause = ""

    # Set up query.
    no_anno = [get_annotation_clause(anno, has_anno=False) for anno in annotations]
    return "%s WHERE (type = '%s') AND (impact_severity = 'MED' or impact_severity = 'HIGH') " \
           "AND (allele_bal >= %f) AND %s AND %s %s" \
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
                #print ("******************* Result already in variants:")
                #print (result.split('\t')[0:3])
            count += 1
        return count

    # Get variants in hotspots.
    for anno in annotations:
        results = get_query_results(gemini_db, get_hotspot_variants(anno), out_format)
        count = add_variants(results)
        var_counts.append(count)

    # Get novel counts.
    for var_type in ['snp', 'indel']:
        results = get_query_results(gemini_db, get_novel_query(annotations, var_type), out_format)
        count = add_variants(results)
        var_counts.append(count)
    return var_counts, variants

def get_query_results(gemini_db, query, out_format):
    """
    Returns results of query.
    """

    gemini = GeminiQuery(gemini_db, out_format=out_format)
    gemini.run(query)
    return gemini


if __name__ == "__main__":
    # Argument setup.
    parser = argparse.ArgumentParser()
    parser.add_argument("--header", action="store_true", help="Print header?")
    parser.add_argument("gemini_db", help="Gemini database to query")
    parser.add_argument("annotations", help="Annotations to query for")
    parser.add_argument("output_vcf", help="Write variants to this file")

    # Argument parsing.
    args = parser.parse_args()
    annotations = args.annotations.split(',')
    gemini_db = args.gemini_db
    db_name, ext = os.path.splitext( os.path.split(gemini_db)[1] )
    output = open(args.output_vcf, 'w')

    # Set up VCF formatter.
    simple_struct = namedtuple('Simple', 'db')
    vcf_format = VCFRowFormat(simple_struct(db=gemini_db))
    print (vcf_format.header(None), file=output)

    # Get counts, variants.
    counts, variants = get_somatic_variants(gemini_db, annotations, out_format=vcf_format)

    # Sort variants by chrom and start position.
    def get_chrom_and_start(v):
        fields = v.split('\t')
        return (int(fields[0][3:]), int(fields[1]))
    variants.sort(key=get_chrom_and_start)
    
    # Print header?
    if args.header:
        print('#Sample %s Novel_SNPS Novel_Indels' % ' '.join(annotations))

    # Print output of variant counts.
    print( '%s %s' % ( db_name, ' '.join(str(c) for c in counts) ) )

    # Print somatic variants VCF.
    for variant in variants:
        print(variant, file=output)

    output.close()

