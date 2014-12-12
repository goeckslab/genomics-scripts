#!/bin/sh
# 
# Create a chromosome.len file from a BAM file header.
# Usage: bam_to_chrom_len.sh <input.bam> > genome.len

# 1. View the header.
# 2. Seleect sequence metadata.
# 3. Split by colon.
# 4. Split by space and print contig name, length in tab-separated performance.
samtools view -H $1 | grep ^@SQ | awk -F ':' '{print $2, $3}' | awk '{printf "%s\t%s\n", $1, $3}'
