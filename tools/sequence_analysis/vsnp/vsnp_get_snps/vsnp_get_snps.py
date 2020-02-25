#!/usr/bin/env python

# Collect quality parsimonious SNPs from vcf files and output alignment files in fasta format.

import argparse
import os
import pandas
import sys
import time
import vcf
from Bio import Phylo
from collections import OrderedDict
from datetime import datetime

# from filter_finder import Filter_Finder

INPUT_VCF_DIR = 'input_vcf_dir'
INPUT_ZC_VCF_DIR = 'input_zc_vcf_dir'
OUTPUT_JSON_SNPS_DIR = 'output_json_snps_dir'
OUTPUT_SNPS_DIR = 'output_snps_dir'


def get_time_stamp():
    return datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H-%M-%S')


class GetSnps:
    def __init__(self, vcf_files, reference, excel_grouper_file, gbk_file, filter_finder,
                 no_filters, all_isolates, ac, mq_val, n_threshold, qual_threshold, output_log):
        self.ac = ac
        self.all_isolates = all_isolates
        # Filter based on the contents of an Excel file.
        self.excel_grouper_file = excel_grouper_file
        self.filter_finder = filter_finder
        # Use Genbank file
        self.gbk_file = gbk_file
        # This will be populated from the columns
        # in the Excel filter file if it is used.
        self.groups = []
        self.mq_val = mq_val
        self.n_threshold = n_threshold
        self.no_filters = no_filters
        # Output process log file handle.
        self.olfh = open(output_log, "w")
        self.qual_threshold = qual_threshold
        # A collection of zero coverage filtered vcf.
        self.vcf_files = vcf_files
        self.olfh.write("Time started: %s\n" % str(get_time_stamp()))
        self.olfh.write("Number VCF inputs: %d\n" % len(self.vcf_files))
        self.olfh.write("Reference: %s\n" % str(reference))
        self.olfh.write("All isolates: %s\n" % str(all_isolates))

    def decide_snps(self, filename):
        # Find the SNPs in a vcf file to produce a pandas data
        # frame and a dictionary containing sample map qualities.
        self.olfh.write("\n%s - Started decide_snps\n" % get_time_stamp())
        sample_map_qualities = {}
        # Eliminate the path.
        file_name_base = os.path.basename(filename)
        # Eliminate the extension.
        file_name_base = os.path.splitext(file_name_base)[0]
        self.olfh.write("\nfile: %s\n" % str(file_name_base))
        vcf_reader = vcf.Reader(open(filename, 'r'))
        sample_dict = {}
        for record in vcf_reader:
            alt = str(record.ALT[0])
            record_position = "%s:%s" % (str(record.CHROM), str(record.POS))
            if record_position in self.all_positions:
                if alt == "None":
                    sample_dict.update({record_position: "-"})
                else:
                    # Not sure this is the best place to capture MQM average
                    # may be faster after parsimony SNPs are decided, but
                    # then it will require opening the files again.
                    # On rare occassions MQM gets called "NaN", thus passing
                    # a string when a number is expected when calculating average.
                    mq_val = self.get_mq_val(record.INFO, filename)
                    if str(mq_val).lower() not in ["nan", "na", "inf"]:
                        sample_map_qualities.update({record_position: mq_val})
                    # Add parameters here to change what each vcf represents.
                    # SNP is represented in table, now how will the vcf represent
                    # the called position alt != "None", which means a deletion
                    # as alt is not record.FILTER, or rather passed.
                    len_alt = len(alt)
                    if len_alt == 1:
                        qual_val = self.val_as_int(record.QUAL)
                        ac = record.INFO['AC'][0]
                        ref = str(record.REF[0])
                        if ac == 2 and qual_val > self.n_threshold:
                            sample_dict.update({record_position: alt})
                        elif ac == 1 and qual_val > self.n_threshold:
                            alt_ref = alt + ref
                            if alt_ref == "AG":
                                sample_dict.update({record_position: "R"})
                            elif alt_ref == "CT":
                                sample_dict.update({record_position: "Y"})
                            elif alt_ref == "GC":
                                sample_dict.update({record_position: "S"})
                            elif alt_ref == "AT":
                                sample_dict.update({record_position: "W"})
                            elif alt_ref == "GT":
                                sample_dict.update({record_position: "K"})
                            elif alt_ref == "AC":
                                sample_dict.update({record_position: "M"})
                            elif alt_ref == "GA":
                                sample_dict.update({record_position: "R"})
                            elif alt_ref == "TC":
                                sample_dict.update({record_position: "Y"})
                            elif alt_ref == "CG":
                                sample_dict.update({record_position: "S"})
                            elif alt_ref == "TA":
                                sample_dict.update({record_position: "W"})
                            elif alt_ref == "TG":
                                sample_dict.update({record_position: "K"})
                            elif alt_ref == "CA":
                                sample_dict.update({record_position: "M"})
                            else:
                                sample_dict.update({record_position: "N"})
                            # Poor calls
                        elif qual_val <= 50:
                            sample_dict.update({record_position: ref})
                        elif qual_val <= self.n_threshold:
                            sample_dict.update({record_position: "N"})
                        else:
                            # Insurance -- Will still report on a possible
                            # SNP even if missed with above statement
                            sample_dict.update({record_position: ref})
        # Merge dictionaries and order
        merge_dict = {}
        # abs_pos:REF
        merge_dict.update(self.all_positions)
        # abs_pos:ALT replacing self.all_positions, because keys must be unique
        merge_dict.update(sample_dict)
        sample_df = pandas.DataFrame(merge_dict, index=[file_name_base])
        return sample_df, file_name_base, sample_map_qualities

    def df_to_fasta(self, parsimonious_df, group, output_fasta):
        # Generate SNP alignment file from the parsimonious_df
        # data frame.  If using an excel filter, group will not
        # be None, but output_fasta will be None.
        self.olfh.write("\n%s - Started df_to_fasta\n" % get_time_stamp())
        if group is None:
            snps_file = output_fasta
        else:
            snps_file = os.path.join(OUTPUT_SNPS_DIR, "%s_parsimonious_snps.fasta" % group)
        test_duplicates = []
        with open(snps_file, 'w') as fh:
            for index, row in parsimonious_df.iterrows():
                test_duplicates.append(row.name)
                if test_duplicates.count(row.name) < 2:
                    fh.write(">%s\n" % row.name)
                    for pos in row:
                        fh.write("%s" % str(pos))
                    fh.write("\n")

    def find_initial_positions(self, filename):
        # Find SNP positions in a vcf file.
        self.olfh.write("\n%s - Started find_initial_positions\n" % get_time_stamp())
        found_positions = {}
        found_positions_mix = {}
        try:
            vcf_reader = vcf.Reader(open(filename, 'r'))
            try:
                for record in vcf_reader:
                    qual_val = self.val_as_int(record.QUAL)
                    chrom = record.CHROM
                    position = record.POS
                    absolute_position = "%s:%s" % (str(chrom), str(position))
                    alt = str(record.ALT[0])
                    if alt != "None":
                        mq_val = self.get_mq_val(record.INFO, filename)
                        ac = record.INFO['AC'][0]
                        len_ref = len(record.REF)
                        if ac == self.ac and len_ref == 1 and qual_val > self.qual_threshold and mq_val > self.mq_val:
                            found_positions.update({absolute_position: record.REF})
                        if ac == 1 and len_ref == 1 and qual_val > self.qual_threshold and mq_val > self.mq_val:
                            found_positions_mix.update({absolute_position: record.REF})
                return filename, found_positions, found_positions_mix
            except (ZeroDivisionError, ValueError, UnboundLocalError, TypeError) as e:
                self.olfh.write("\nException thrown parsing record in file %s: %s\n" % (filename, str(e)))
                return filename, {'': ''}, {'': ''}
        except (SyntaxError, AttributeError) as e:
            self.olfh.write("\nException thrown by vcf.Reader attempting to read file %s: %s\n" % (filename, str(e)))
            return filename, {'': ''}, {'': ''}

    def gather_and_filter(self, prefilter_df, group, output_fasta, output_json_snps):
        # Group and filter adata frame of SNPs.  If an
        # excel file is not used, group will be None and
        # output_fasta and output_json_snps will not be
        # None, and vice versa.
        self.olfh.write("\n%s - Started gather_and_filter\n" % get_time_stamp())
        if self.excel_grouper_file is None or self.no_filters:
            filtered_all_df = prefilter_df
            sheet_names = None
        elif self.excel_grouper_file:
            # The value of group is not None.
            # Filter positions to be removed from all.
            xl = pandas.ExcelFile(self.excel_grouper_file)
            sheet_names = xl.sheet_names
            # Use the first column to filter "all" postions.
            exclusion_list_all = self.get_position_list(sheet_names, 0)
            exclusion_list_group = self.get_position_list(sheet_names, group)
            exclusion_list = exclusion_list_all + exclusion_list_group
            # Filters for all applied.
            filtered_all_df = prefilter_df.drop(columns=exclusion_list, errors='ignore')
        parsimonious_df = self.get_parsimonious_df(filtered_all_df)
        if group is not None:
            json_snps_file = os.path.join(OUTPUT_JSON_SNPS_DIR, "%s_snps.json" % group)
        else:
            json_snps_file = output_json_snps
        parsimonious_df.to_json(json_snps_file, orient='split')
        self.df_to_fasta(parsimonious_df, group, output_fasta)
        samples_number, columns = parsimonious_df.shape
        if samples_number < 4:
            msg = "Too few samples to build tree"
            if group is not None:
                msg = "%s for group: %s" % (msg, group)
            self.olfh.write("%s\n" % msg)

    def get_mq_val(self, record_info, filename):
        # Get the MQ (gatk) or MQM (freebayes) value
        # from the record.INFO component of the vcf file.
        try:
            mq_val = record_info['MQM']
            return self.return_val(mq_val)
        except Exception:
            try:
                mq_val = record_info['MQ']
                return self.return_val(mq_val)
            except Exception:
                msg = "Invalid or unsupported vcf header %s in file: %s\n" % (str(record_info), filename)
                sys.exit(msg)

    def get_parsimonious_df(self, filtered_all_df):
        # Get the parsimonious SNPs data frame
        # from a data frame of filtered SNPs.
        self.olfh.write("\n%s - Started get_parsimonious_df\n" % get_time_stamp())
        ref_series = filtered_all_df.loc['root']
        # In all_vcf root needs to be removed.
        filtered_all_df = filtered_all_df.drop(['root'])
        parsimony = filtered_all_df.loc[:, (filtered_all_df != filtered_all_df.iloc[0]).any()]
        parsimony_positions = list(parsimony)
        parse_df = filtered_all_df[parsimony_positions]
        ref_df = ref_series.to_frame()
        ref_df = ref_df.T
        parsimonious_df = pandas.concat([parse_df, ref_df], join='inner')
        return parsimonious_df

    def get_position_list(self, sheet_names, group):
        # Get a list of positions defined by an excel file.
        self.olfh.write("\n%s - Started get_position_list\n" % get_time_stamp())
        exclusion_list = []
        try:
            filter_to_all = pandas.read_excel(self.excel_grouper_file, header=1, usecols=[group])
            for value in filter_to_all.values:
                value = str(value[0])
                if "-" not in value.split(":")[-1]:
                    exclusion_list.append(value)
                elif "-" in value:
                    try:
                        chrom, sequence_range = value.split(":")
                    except Exception as e:
                        sys.exit(str(e))
                    value = sequence_range.split("-")
                    for position in range(int(value[0].replace(',', '')), int(value[1].replace(',', '')) + 1):
                        exclusion_list.append(chrom + ":" + str(position))
            return exclusion_list
        except ValueError:
            exclusion_list = []
            return exclusion_list

    def val_as_int(self, val):
        # Handle integer value conversion.
        try:
            return int(val)
        except TypeError:
            # val is likely None here.
            return 0

    def get_snps(self, output_json_avg_mq, output_json_snps=None, group=None, output_fasta=None):
        # Parse all vcf files to accumulate SNPs into a
        # data frame.  If group is None, output_fasta will
        # not be None and vice versa.
        self.olfh.write("\n%s - Started get_snps\n" % get_time_stamp())
        if self.filter_finder:
            # Process Excel filter file.
            self.self.olfh.write(f'\nFinding filters...\n')
            # TODO: fix this...
            # filter_finder = Filter_Finder(self.excel_grouper_file)
            # filter_finder.filter_finder()
        all_positions = {}
        df_list = []
        self.olfh.write("self.vcf_files: %s\n" % str(self.vcf_files))
        for filename in self.vcf_files:
            filename2, found_positions, found_positions_mix = self.find_initial_positions(filename)
            all_positions.update(found_positions)
        # Order before adding to file to match
        # with ordering of individual samples.
        # all_positions is abs_pos:REF
        self.all_positions = OrderedDict(sorted(all_positions.items()))
        ref_positions_df = pandas.DataFrame(self.all_positions, index=['root'])
        all_map_qualities = {}
        for filename in self.vcf_files:
            sample_df, file_name_base, sample_map_qualities = self.decide_snps(filename)
            df_list.append(sample_df)
            all_map_qualities.update({file_name_base: sample_map_qualities})
        all_sample_df = pandas.concat(df_list)
        # All positions have now been selected for each sample,
        # so select parisomony informative SNPs.  This removes
        # columns where all fields are the same.
        # Add reference to top row.
        prefilter_df = pandas.concat([ref_positions_df, all_sample_df], join='inner')
        all_mq_df = pandas.DataFrame.from_dict(all_map_qualities)
        mq_averages = all_mq_df.mean(axis=1).astype(int)
        mq_averages.to_json(output_json_avg_mq, orient='split')
        self.gather_and_filter(prefilter_df, group, output_fasta, output_json_snps)

    def group_vcfs(self):
        # Parse an excel file to produce a
        # grouping dictionary for filtering SNPs.
        self.olfh.write("\n%s - Started group_vcfs\n" % get_time_stamp())
        group_names = []
        xl = pandas.ExcelFile(self.excel_grouper_file)
        sheet_names = xl.sheet_names
        ws = pandas.read_excel(self.excel_grouper_file, sheet_name=sheet_names[0])
        defining_snps = ws.iloc[0]
        defsnp_iterator = iter(defining_snps.iteritems())
        next(defsnp_iterator)
        defining_snps = {}
        for abs_pos, group in defsnp_iterator:
            group_names.append(group)
            # Make defining snp/group dict.
            defining_snps[abs_pos] = group
        samples_groups_dict = {}
        for vcf_file in self.vcf_files:
            filename, found_positions, found_positions_mix = self.find_initial_positions(vcf_file)
            filename = os.path.basename(filename)
            groups = []
            try:
                for abs_position in list(defining_snps.keys() & (found_positions.keys() | found_positions_mix.keys())):
                    group = defining_snps[abs_position]
                    groups.append(group)
                if groups:
                    groups = sorted(groups)
                    samples_groups_dict[vcf_file] = groups
                    self.groups = groups
                else:
                    samples_groups_dict[filename] = []
            except TypeError as e:
                self.olfh.write("Exception processing file %s when grouping vcfs: %s\n" % (filename, str(e)))
                samples_groups_dict[filename] = []
        return samples_groups_dict

    def return_val(self, val):
        # Handle element and single-element list values.
        if isinstance(val, list):
            return val[0]
        return val

    def tree_to_svg(self, tree_file, svg_file):
        # Convert a phylogenetic tree to an svg output.
        self.olfh.write("\n%s - Started tree_to_svg\n" % get_time_stamp())

        def get_label(leaf):
            return leaf.name

        tree = Phylo.read(tree_file, 'newick')
        tree.ladderize()
        Phylo.draw(tree, label_func=get_label, do_show=False)
        # pylab.axis('off')
        # pylab.savefig(svg_file, format='svg', bbox_inches='tight', dpi=500)


parser = argparse.ArgumentParser()

parser.add_argument('--all_isolates', action='store_true', dest='all_isolates', required=False, default=False, help='Create table with all isolates'),
parser.add_argument('--filter_finder', action='store_true', dest='filter_finder', required=False, default=False, help='Write possible positions to filter to text file'),
parser.add_argument('--excel_grouper_file', action='store', dest='excel_grouper_file', required=False, default=None, help='Optional Excel filter file'),
parser.add_argument('--gbk_file', action='store', dest='gbk_file', required=False, default=None, help='Optional gbk file'),
parser.add_argument('--no_filters', action='store_true', dest='no_filters', default=False, help='Run without applying filters'),
parser.add_argument('--output_fasta', action='store', dest='output_fasta', required=False, default=None, help='Single output SNPs alignment fasta file if not Excel filtering'),
parser.add_argument('--output_json_avg_mq', action='store', dest='output_json_avg_mq', help='Single output average mq json file if not Excel filtering'),
parser.add_argument('--output_json_snps', action='store', dest='output_json_snps', required=False, default=None, help='Single output parsimonious SNPs json file if not Excel filtering'),
parser.add_argument('--output_log', action='store', dest='output_log', help='Output log file'),
parser.add_argument('--reference', action='store', dest='reference', help='Reference file'),
parser.add_argument('--subset', action='store_true', dest='subset', required=False, default=False, help='Create trees with a subset of sample that represent the whole'),

args = parser.parse_args()

# Initializations - TODO: should these be passed in as command line args?
ac = 2
mq_val = 56
n_threshold = 50
qual_threshold = 150
# Build the list of vcf files against which
# the current run will be analyzed.  This set
# of files will be in a directory named input_vcf_dir.
vcf_files = []
for file_name in os.listdir(INPUT_VCF_DIR):
    file_path = os.path.abspath(os.path.join(INPUT_VCF_DIR, file_name))
    vcf_files.append(file_path)
# Add the zero coverage vcf files from the current run.
for file_name in os.listdir(INPUT_ZC_VCF_DIR):
    file_path = os.path.abspath(os.path.join(INPUT_ZC_VCF_DIR, file_name))
    vcf_files.append(file_path)
snp_finder = GetSnps(vcf_files, args.reference, args.excel_grouper_file, args.gbk_file, args.filter_finder, args.no_filters,
                     args.all_isolates, ac, mq_val, n_threshold, qual_threshold, args.output_log)

if args.excel_grouper_file is not None:
    # Parse the Excel file to detemine groups for filtering.
    samples_groups_dict = snp_finder.group_vcfs()
    for group in snp_finder.groups:
        snp_finder.get_snps(args.output_json_avg_mq, output_json_snps=None, group=group, output_fasta=None)
if args.all_isolates or args.subset or args.excel_grouper_file is None:
    snp_finder.get_snps(args.output_json_avg_mq,
                        output_json_snps=args.output_json_snps,
                        group=None,
                        output_fasta=args.output_fasta)
snp_finder.olfh.write("\nTime finished: %s\n\n" % get_time_stamp())
snp_finder.olfh.close()
