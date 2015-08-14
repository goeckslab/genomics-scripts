#
# Script for simple statistical testing of replicates vs. non-replicates.
#


import sys
import re
import itertools
import pandas
from scipy import stats
from numpy import mean, var, std, histogram

class Sample:
    def __init__(self, name=None, num_vars=None):
        self.name = name
        self.num_vars = num_vars


def get_concordances(data, reduce_fn=mean, weighted=False):
    """
    Reads variant genotype data and returns a tuple of <replicate_concordance> <nonreplicate_concordance> where each
    element has concordances of replicates and non-replicates.

    reduce_fn - reduce function used to reduce % concordance of two files into a single number
    weighted - if True, weight concordance by number of shared variants; if False, use only % concordance

    For now, data format has a header line and corresponding data and is tab separated:
    variant_id chrom start gt1 gt2 ... gtN
    12 chr1 1234567 ./. A/T ... A/T
    ...
    """

    # TODO: should be able to provide the data programmatically from GEMINI
    """
    query = "SELECT variant_id, chrom, start, (gts).(*) FROM variants"
    structured_variants = []
    results = gem_ops.get_query_results( gemini_db, query) )
    # Add header to results and feed in as data below.
    """

    genotypes_df = pandas.DataFrame.from_csv(data, sep='\t')
    sample_names = list(genotypes_df.columns.values)[2:]
    sample_names.sort()

    # Get number of variants for each sample and create Sample objects.
    samples = []
    for name in sample_names:
        num_vars = len(genotypes_df.loc[genotypes_df[name] != './.', name])
        samples.append( Sample(name, num_vars) )


    # Compute concordance amongst samples.
    replicate_concordance = []
    nonreplicate_concordance = []
    for sample1, sample2 in itertools.combinations(samples, 2):
        # Get shared variants between samples and overall concordance.
        num_shared_vars = len ( genotypes_df.loc[(genotypes_df[sample1.name] != './.') & \
                                                 (genotypes_df[sample1.name] == genotypes_df[sample2.name]), \
                                                 sample2.name] )
        pct1 = num_shared_vars/float(sample1.num_vars) * 100
        pct2 = num_shared_vars/float(sample2.num_vars) * 100
        concordance = reduce_fn([pct1, pct2])

        #print sample1.name, sample2.name, num_shared_vars, pct1, pct2

        # If weighting, weight concordance by number of concordant variants; multiple by 0.01
        # to convert concordance % to # of variants.
        if weighted:
            concordance *= num_shared_vars * 0.01

        # If samples are replicates, add to replicates; otherwise add to
        # nonreplicates.
        if sample1.name.replace('FirstBatch', '') == sample2.name.replace('Repeats', ''):
            # Replicate.
            #print sample1.name, sample2.name, num_shared_vars, pct1, pct2
            replicate_concordance.append(concordance)
        else:
            # Not replicates.
            nonreplicate_concordance.append(concordance)

    return replicate_concordance, nonreplicate_concordance


if __name__ == "__main__":

    # Read data.
    replicate_concordance, nonreplicate_concordance = get_concordances(sys.stdin, weighted=False)

    # Run stats. NOTE: neither set of concordances is normally distributed, but the CLT
    # states that with a large enough sample, t-tests can be used. Alternatively,
    # Mann-Whitney will work as well.
    #
    # References:
    # http://stats.stackexchange.com/questions/15664/how-to-test-for-differences-between-two-group-means-when-the-data-is-not-normall
    # http://stats.stackexchange.com/questions/71452/testing-for-significance-between-means-having-one-normal-distributed-sample-and
    # http://stats.stackexchange.com/questions/9573/t-test-for-non-normal-when-n50

    # Print descriptive statistics.
    print "Counts (rep, non-rep): %i, %i" % ( len(replicate_concordance), len(nonreplicate_concordance) )
    print "Mins (rep, non-rep): %.3f, %.3f" % ( min(replicate_concordance), min(nonreplicate_concordance) )
    print "Maxs (rep, non-rep): %.3f, %.3f" % ( max(replicate_concordance), max(nonreplicate_concordance) )
    print "Means (rep, non-rep): %.3f, %.3f" % ( mean(replicate_concordance), mean(nonreplicate_concordance) )
    print "Variance (rep, non-rep): %.3f, %.3f" % ( var(replicate_concordance), var(nonreplicate_concordance) )
    print "Replicates histogram: ", histogram(replicate_concordance)
    print "Nonreplicates histogram: ", histogram(nonreplicate_concordance)

    # Unpaired t-test assuming equal variances.
    two_sample = stats.ttest_ind(replicate_concordance, nonreplicate_concordance)
    print "The t-statistic is %.3f and the p-value is %.3E." % two_sample

    # Unpaired t-test assuming unequal population variances.
    two_sample_diff_var = stats.ttest_ind(replicate_concordance, nonreplicate_concordance, equal_var=False)
    print "If we assume unequal variances than the t-statistic is %.3f and the p-value is %.3E." % two_sample_diff_var

    # Mann Whitney U test.
    mann_whitney = stats.mannwhitneyu(replicate_concordance, nonreplicate_concordance)
    print "Mann Whitney U test: MW-statistic is %.3f and p-value is %.3E" % mann_whitney
