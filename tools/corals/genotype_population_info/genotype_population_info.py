#!/usr/bin/env python
"""
Generate the genotype_population_info.txt file by parsing the information
from a Affymetrix 96 well plate CSV file and an associated VCF file.
"""
import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument('--input_csv', dest='input_csv', help='Affymetrix 96 well plate file')
parser.add_argument('--input_vcf', dest='input_vcf', help='Input VCF file')
parser.add_argument('--output', dest='output', help='Output dataset'),
args = parser.parse_args()

# Parse the input_vcf file, looking for the first line
# that starts with the string "#CHROM"
with open(args.input_vcf, "r") as vcfh:
    for line in vcfh:
        if not line.startswith("#CHROM"):
            continue
        line = line.rstrip("\r\n")
        # Example line:
        # #CHROM  13704   13706   13708   13736   13748   13762   13782
        items = line.split("\t")
        sample_list = items[8:]
        break

# Parse the input_csv file to get the region for for
# each sample_id in the sample_list.  Initialize the
# region_list to be the same as the sample_list to ensure
# the same length.
region_list = [x for x in sample_list]
with open(args.input_csv, "r") as csvh:
    for i, line in enumerate(csvh):
        if i == 0:
            # Skip the header.
            continue
        line = line.rstrip('\r\n')
        items = line.split(',')
        csv_sample_id = items[0]
        csv_region = items[9]
        # Make sure the csv_sample_id is in the sample_list.
        try:
            loc = sample_list.index(csv_sample_id)
            region_list[loc] = csv_region
        except Exception:
            pass

# The output file will consist of columns:
# Item #, Sample ID, Region
with open(args.output, "w") as outfh:
    for i, sample_id in enumerate(sample_list):
        outfh.write("%d\t%s\t%s\n" % (i, sample_id, region_list[1]))
