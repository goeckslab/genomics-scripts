import argparse
import os
from gemini import GeminiQuery


COMMON_DATABASES = ["1kg", "exac", "esp"]
ANNOTATIONS = ["TCGA_LUAD"]

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
    return ' AND '.join( ['(%s)' % get_in_and_aff_clause(db, 0.01) for db in COMMON_DATABASES ] )

def get_annotation_and_no_common_clause(anno_col_name, aaf=0.01):
    """
    Returns clause for selecting variants with an annotation but not appearing in common databases at a given alternate allele frequency.
    """
    return '%s=1 AND %s' % (anno_col_name, get_no_common_vars_clause())

def get_count(gemini, query):
    gemini.run(query)
    for i, row in enumerate(gemini):
        pass
    return i + 1

def somatic_query(sample_db):
    gq = GeminiQuery(sample_db)
    query = "select chrom, start, end, ref, alt, aaf_1kg_all, aaf_exac_all, aaf_esp_all, TCGA_LUAD from variants WHERE " + get_annotation_and_no_common_clause("TCGA_LUAD")
    print query
    gq.run(query)
    for row in gq:
        print row
    return get_count(GeminiQuery(sample_db), query)

if __name__ == "__main__":
    # Argument parsing.
    parser = argparse.ArgumentParser()
    parser.add_argument("sample_db", help="Gemini database to query")
    args = parser.parse_args()

    # Do queries.
    db_name, ext = os.path.splitext( os.path.split(args.sample_db)[1] )
    counts = [somatic_query(args.sample_db)]

    # Print output.
    print '%s %s' % ( db_name, ' '.join(str(c) for c in counts) ) 


