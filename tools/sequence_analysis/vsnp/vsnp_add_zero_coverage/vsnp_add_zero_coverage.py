#!/usr/bin/env python

import argparse
import shutil
import re
import os
import pandas
import pysam
from numpy import mean
from Bio import SeqIO

HEADERFILE = "v_header.csv"
BODYFILE = "v_annotated_body.csv"

parser = argparse.ArgumentParser()

parser.add_argument('--bam_input', action='store', dest='bam_input', help='Input BAM file')
parser.add_argument('--vcf_input', action='store', dest='vcf_input', help='Input VCF file with poor positions removed')
parser.add_argument('--output_metrics', action='store', dest='output_metrics', help='Output metrics text file')
parser.add_argument('--output_vcf', action='store', dest='output_vcf', help='Output VCF file')
parser.add_argument('--reference', action='store', dest='reference', help='Reference dataset')

args = parser.parse_args()

# Create a coverage dictionary.
coverage_dict = {}
coverage_list = pysam.depth(args.bam_input, split_lines=True)
for line in coverage_list:
    chrom, position, depth = line.split('\t')
    coverage_dict["%s-%s" % (chrom, position)] = depth
# Convert it to a data frame.
coverage_df = pandas.DataFrame.from_dict(coverage_dict, orient='index', columns=["depth"])
# Create a zero coverage dictionary.
zero_dict = {}
for record in SeqIO.parse(args.reference, "fasta"):
    chrom = record.id
    total_len = len(record.seq)
    for pos in list(range(1, total_len + 1)):
        zero_dict["%s-%s" % (str(chrom), str(pos))] = 0
# Convert it to a data frame with depth_x
# and depth_y columns - index is NaN.
zero_df = pandas.DataFrame.from_dict(zero_dict, orient='index', columns=["depth"])
coverage_df = zero_df.merge(coverage_df, left_index=True, right_index=True, how='outer')
# depth_x "0" column no longer needed.
coverage_df = coverage_df.drop(columns=['depth_x'])
coverage_df = coverage_df.rename(columns={'depth_y': 'depth'})
# Covert the NaN to 0 coverage and get some metrics.
coverage_df = coverage_df.fillna(0)
coverage_df['depth'] = coverage_df['depth'].apply(int)
total_length = len(coverage_df)
average_coverage = coverage_df['depth'].mean()
zero_df = coverage_df[coverage_df['depth'] == 0]
total_zero_coverage = len(zero_df)
total_coverage = total_length - total_zero_coverage
genome_coverage = "{:.2%}".format(total_coverage / total_length)
# Read the input.
column_names = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "Sample"]
vcf_df = pandas.read_csv(args.vcf_input, sep='\t', header=None, names=column_names, comment='#')
good_snp_count = len(vcf_df[(vcf_df['ALT'].str.len() == 1) & (vcf_df['REF'].str.len() == 1) & (vcf_df['QUAL'] > 150)])
if total_zero_coverage > 0:
    header_file = HEADERFILE
    with open(header_file, 'w') as outfile:
        with open(args.vcf_input) as infile:
            for line in infile:
                if re.search('^#', line):
                    outfile.write("%s" % line)
    vcf_df_snp = vcf_df[vcf_df['REF'].str.len() == 1]
    vcf_df_snp = vcf_df_snp[vcf_df_snp['ALT'].str.len() == 1]
    vcf_df_snp['ABS_VALUE'] = vcf_df_snp['CHROM'].map(str) + "-" + vcf_df_snp['POS'].map(str)
    vcf_df_snp = vcf_df_snp.set_index('ABS_VALUE')
    cat_df = pandas.concat([vcf_df_snp, zero_df], axis=1, sort=False)
    cat_df = cat_df.drop(columns=['CHROM', 'POS', 'depth'])
    cat_df[['ID', 'ALT', 'QUAL', 'FILTER', 'INFO']] = cat_df[['ID', 'ALT', 'QUAL', 'FILTER', 'INFO']].fillna('.')
    cat_df['REF'] = cat_df['REF'].fillna('N')
    cat_df['FORMAT'] = cat_df['FORMAT'].fillna('GT')
    cat_df['Sample'] = cat_df['Sample'].fillna('./.')
    cat_df['temp'] = cat_df.index.str.rsplit('-', n=1)
    cat_df[['CHROM', 'POS']] = pandas.DataFrame(cat_df.temp.values.tolist(), index=cat_df.index)
    cat_df = cat_df[['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'Sample']]
    cat_df['POS'] = cat_df['POS'].astype(int)
    cat_df = cat_df.sort_values(['CHROM', 'POS'])
    body_file = BODYFILE
    cat_df.to_csv(body_file, sep='\t', header=False, index=False)
    with open(args.output_vcf, "w") as outfile:
        for cf in [HEADERFILE, BODYFILE]:
            with open(cf, "r") as infile:
                for line in infile:
                    outfile.write("%s" % line)
else:
    shutil.copyfile(args.vcf_input, args.output_vcf)

metrics_columns = ["File", "Number of Good SNPs", "Average Coverage", "Genome Coverage"]
bam_metrics = [os.path.basename(args.bam_input), "", "%4f" % average_coverage, genome_coverage]
vcf_metrics = [os.path.basename(args.vcf_input), str(good_snp_count), "", ""]
with open(args.output_metrics, "w") as fh:
    fh.write("# %s\n" % "\t".join(metrics_columns))
    fh.write("%s\n" % "\t".join(bam_metrics))
    fh.write("%s\n" % "\t".join(vcf_metrics))
