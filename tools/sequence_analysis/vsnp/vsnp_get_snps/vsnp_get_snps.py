#!/usr/bin/env python

# Collect quality parsimonious SNPs from vcf files and output alignment files in fasta format.

import argparse
import glob
import os
import pandas
import shutil
import sys
import time
import vcf
from collections import OrderedDict
from datetime import datetime

ALL_VCFS_DIR = 'all_vcf'
INPUT_VCF_DIR = 'input_vcf_dir'
OUTPUT_JSON_SNPS_DIR = 'output_json_snps_dir'
OUTPUT_SNPS_DIR = 'output_snps_dir'


def get_time_stamp():
    return datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H-%M-%S')


def setup_all_vcfs(vcf_files, vcf_dirs):
    # Create the all_vcfs directory and link
    # all input vcf files into it for processing.
    os.makedirs(ALL_VCFS_DIR)
    vcf_dirs.append(ALL_VCFS_DIR)
    for vcf_file in vcf_files:
        file_name_base = os.path.basename(vcf_file)
        dst_file = os.path.join(ALL_VCFS_DIR, file_name_base)
        os.symlink(vcf_file, dst_file)
    return vcf_dirs


class GetSnps:

    def __init__(self, vcf_files, reference, excel_grouper_file, gbk_file,
                 all_isolates, ac, mq_val, n_threshold, qual_threshold, output_summary):
        self.ac = ac
        self.all_isolates = all_isolates
        self.all_positions = None
        # Filter based on the contents of an Excel file.
        self.excel_grouper_file = excel_grouper_file
        # Use Genbank file
        self.gbk_file = gbk_file
        self.groups = []
        # This will be populated from the columns
        # in the Excel filter file if it is used.
        self.mq_val = mq_val
        self.n_threshold = n_threshold
        self.olfh = None
        self.qual_threshold = qual_threshold
        self.reference = reference
        self.start_time = get_time_stamp()
        self.timer_start = datetime.now()
        self.vcf_files = vcf_files
        self.initiate_summary(output_summary)

    def bin_input_files(self, filename, samples_groups_dict, defining_snps, inverted_defining_snps, found_positions, found_positions_mix):
        sample_groups_list = []
        table_name = self.get_base_file_name(filename)
        try:
            defining_snp = False
            # Absolute positions in set union of two lists.
            for abs_position in list(defining_snps.keys() & (found_positions.keys() | found_positions_mix.keys())):
                group = defining_snps[abs_position]
                sample_groups_list.append(group)
                self.check_add_group(group)
                if len(list(defining_snps.keys() & found_positions_mix.keys())) > 0:
                    table_name = self.get_base_file_name(filename)
                    table_name = '%s<font color="red">[[MIXED]]</font>' % table_name
                self.copy_file(filename, group)
                defining_snp = True
            if not set(inverted_defining_snps.keys()).intersection(found_positions.keys() | found_positions_mix.keys()):
                for abs_position in list(inverted_defining_snps.keys()):
                    group = inverted_defining_snps[abs_position]
                    sample_groups_list.append(group)
                    self.check_add_group(group)
                    self.copy_file(filename, group)
                    defining_snp = True
            if defining_snp:
                samples_groups_dict[table_name] = sorted(sample_groups_list)
            else:
                samples_groups_dict[table_name] = ['<font color="red">No defining SNP</font>']
        except TypeError as e:
            msg = "<br/>Error processing file %s to generate samples_groups_dict: %s<br/>" % (filename, str(e))
            self.olfh.write(msg)
            samples_groups_dict[table_name] = [msg]
        return samples_groups_dict

    def check_add_group(self, group):
        if group not in self.groups:
            self.groups.append(group)

    def copy_file(self, filename, dir):
        if not os.path.exists(dir):
            os.makedirs(dir)
        shutil.copy(filename, dir)

    def decide_snps(self, filename):
        # Find the SNPs in a vcf file to produce a pandas data
        # frame and a dictionary containing sample map qualities.
        sample_map_qualities = {}
        # Eliminate the path.
        file_name_base = self.get_base_file_name(filename)
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
        # abs_pos:ALT replacing all_positions, because keys must be unique
        merge_dict.update(sample_dict)
        sample_df = pandas.DataFrame(merge_dict, index=[file_name_base])
        return sample_df, file_name_base, sample_map_qualities

    def df_to_fasta(self, parsimonious_df, group):
        # Generate SNP alignment file from the parsimonious_df
        # data frame.
        snps_file = os.path.join(OUTPUT_SNPS_DIR, "%s.fasta" % group)
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
                return found_positions, found_positions_mix
            except (ZeroDivisionError, ValueError, UnboundLocalError, TypeError) as e:
                self.olfh.write("<br/>Error parsing record in file %s: %s<br/>" % (filename, str(e)))
                return {'': ''}, {'': ''}
        except (SyntaxError, AttributeError) as e:
            self.olfh.write("<br/>Error attempting to read file %s: %s<br/>" % (filename, str(e)))
            return {'': ''}, {'': ''}

    def gather_and_filter(self, prefilter_df, group_dir):
        # Group a data frame of SNPs.
        if self.excel_grouper_file is None:
            filtered_all_df = prefilter_df
            sheet_names = None
        else:
            # Filter positions to be removed from all.
            xl = pandas.ExcelFile(self.excel_grouper_file)
            sheet_names = xl.sheet_names
            # Use the first column to filter "all" postions.
            exclusion_list_all = self.get_position_list(sheet_names, 0)
            exclusion_list_group = self.get_position_list(sheet_names, group_dir)
            exclusion_list = exclusion_list_all + exclusion_list_group
            # Filters for all applied.
            filtered_all_df = prefilter_df.drop(columns=exclusion_list, errors='ignore')
        json_snps_file = os.path.join(OUTPUT_JSON_SNPS_DIR, "%s.json" % group_dir)
        parsimonious_df = self.get_parsimonious_df(filtered_all_df)
        samples_number, columns = parsimonious_df.shape
        if samples_number >= 4:
            self.df_to_fasta(parsimonious_df, group_dir)
            parsimonious_df.to_json(json_snps_file, orient='split')
        else:
            msg = "<br/>Too few samples to build tree"
            if group_dir is not None:
                msg = "%s for group: %s" % (msg, group_dir)
            self.olfh.write("%s<br/>\n" % msg)

    def get_base_file_name(self, file_path):
        file_name_base = os.path.basename(file_path)
        # Eliminate the extension.
        return os.path.splitext(file_name_base)[0]

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
        try:
            ref_series = filtered_all_df.loc['root']
            # In all_vcf root needs to be removed.
            filtered_all_df = filtered_all_df.drop(['root'])
        except KeyError:
            pass
        parsimony = filtered_all_df.loc[:, (filtered_all_df != filtered_all_df.iloc[0]).any()]
        parsimony_positions = list(parsimony)
        parse_df = filtered_all_df[parsimony_positions]
        ref_df = ref_series.to_frame()
        ref_df = ref_df.T
        parsimonious_df = pandas.concat([parse_df, ref_df], join='inner')
        return parsimonious_df

    def get_position_list(self, sheet_names, group):
        # Get a list of positions defined by an excel file.
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

    def get_snps(self, group_dir, output_json_avg_mq):
        # Parse all vcf files to accumulate SNPs into a
        # data frame.
        all_positions = {}
        df_list = []
        vcf_files = glob.glob(f'{group_dir}/*')
        for vcf_file in vcf_files:
            try:
                found_positions, found_positions_mix = self.find_initial_positions(vcf_file)
                all_positions.update(found_positions)
            except Exception as e:
                self.olfh.wirte("Error updating the all_positions dictionary when processing file %s:\n%s\n" % (vcf_file, str(e)))
        # Order before adding to file to match
        # with ordering of individual samples.
        # all_positions is abs_pos:REF
        self.all_positions = OrderedDict(sorted(all_positions.items()))
        ref_positions_df = pandas.DataFrame(self.all_positions, index=['root'])
        all_map_qualities = {}
        for vcf_file in vcf_files:
            sample_df, file_name_base, sample_map_qualities = self.decide_snps(vcf_file)
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
        self.gather_and_filter(prefilter_df, group_dir)

    def group_vcfs(self, vcf_files):
        # Parse an excel file to produce a
        # grouping dictionary for filtering SNPs.
        xl = pandas.ExcelFile(self.excel_grouper_file)
        sheet_names = xl.sheet_names
        ws = pandas.read_excel(self.excel_grouper_file, sheet_name=sheet_names[0])
        defining_snps = ws.iloc[0]
        defsnp_iterator = iter(defining_snps.iteritems())
        next(defsnp_iterator)
        defining_snps = {}
        inverted_defining_snps = {}
        for abs_pos, group in defsnp_iterator:
            if '!' in abs_pos:
                inverted_defining_snps[abs_pos.replace('!', '')] = group
            else:
                defining_snps[abs_pos] = group
        samples_groups_dict = {}
        for vcf_file in vcf_files:
            found_positions, found_positions_mix = self.find_initial_positions(vcf_file)
            samples_groups_dict = self.bin_input_files(vcf_file, samples_groups_dict, defining_snps, inverted_defining_snps, found_positions, found_positions_mix)
        # Output summary grouping table.
        self.olfh.write('<br/>')
        self.olfh.write('<b>Groupings with %d listed:</b><br/>\n' % len(samples_groups_dict))
        self.olfh.write('<table  cellpadding="5" cellspaging="5" border="1">\n')
        for key, value in samples_groups_dict.items():
            self.olfh.write('<tr align="left"><th>Sample Name</th>\n')
            self.olfh.write('<td>%s</td>' % key)
            for group in value:
                self.olfh.write('<td>%s</td>\n' % group)
            self.olfh.write('</tr>\n')
        self.olfh.write('</table><br/>\n')

    def initiate_summary(self, output_summary):
        # Output summary file handle.
        self.olfh = open(output_summary, "w")
        self.olfh.write('<html>\n')
        self.olfh.write('<head></head>\n')
        self.olfh.write('<body style=\"font-size:12px;">')
        self.olfh.write("<b>Time started:</b> %s<br/>" % str(get_time_stamp()))
        self.olfh.write("<b>Number of VCF inputs:</b> %d<br/>" % len(self.vcf_files))
        self.olfh.write("<b>Reference:</b> %s<br/>" % str(self.reference))
        self.olfh.write("<b>All isolates:</b> %s<br/>" % str(self.all_isolates))

    def return_val(self, val, index=0):
        # Handle element and single-element list values.
        if isinstance(val, list):
            return val[index]
        return val

    def val_as_int(self, val):
        # Handle integer value conversion.
        try:
            return int(val)
        except TypeError:
            # val is likely None here.
            return 0


parser = argparse.ArgumentParser()

parser.add_argument('--all_isolates', action='store', dest='all_isolates', required=False, default="No", help='Create table with all isolates'),
parser.add_argument('--excel_grouper_file', action='store', dest='excel_grouper_file', required=False, default=None, help='Optional Excel filter file'),
parser.add_argument('--gbk_file', action='store', dest='gbk_file', required=False, default=None, help='Optional gbk file'),
parser.add_argument('--output_json_avg_mq', action='store', dest='output_json_avg_mq', help='Single output average mq json file if not Excel filtering'),
parser.add_argument('--output_summary', action='store', dest='output_summary', help='Output summary html file'),
parser.add_argument('--reference', action='store', dest='reference', help='Reference file'),

args = parser.parse_args()

# Initializations - TODO: should these be passed in as command line args?
ac = 2
mq_val = 56
n_threshold = 50
qual_threshold = 150
# Build the list of sample vcf files for the current run.
vcf_files = []
for file_name in os.listdir(INPUT_VCF_DIR):
    file_path = os.path.abspath(os.path.join(INPUT_VCF_DIR, file_name))
    vcf_files.append(file_path)
# Initialize the snp_finder object.
snp_finder = GetSnps(vcf_files, args.reference, args.excel_grouper_file, args.gbk_file,
                     args.all_isolates, ac, mq_val, n_threshold, qual_threshold, args.output_summary)
# Initialize the set of directories containiing vcf files for analysis.
vcf_dirs = []
if args.excel_grouper_file is None:
    vcf_dirs = setup_all_vcfs(vcf_files, vcf_dirs)
    for vcf_dir in vcf_dirs:
        snp_finder.get_snps(vcf_dir, args.output_json_avg_mq)
else:
    if args.all_isolates.lower() == "yes":
        vcf_dirs = setup_all_vcfs(vcf_files, vcf_dirs)
    # Parse the Excel file to detemine groups for filtering.
    snp_finder.group_vcfs(vcf_files)
    # Append the list of group directories created by
    # the above call to the set of directories containing
    # vcf files for analysis
    group_dirs = [d for d in os.listdir(os.getcwd()) if os.path.isdir(d) and d in snp_finder.groups]
    vcf_dirs.extend(group_dirs)
    for vcf_dir in vcf_dirs:
        snp_finder.get_snps(vcf_dir, args.output_json_avg_mq)
# Finish summary log.
snp_finder.olfh.write("<br/><b>Time finished:</b> %s<br/>\n" % get_time_stamp())
total_run_time = datetime.now() - snp_finder.timer_start
snp_finder.olfh.write("<br/><b>Total run time:</b> %s<br/>\n" % str(total_run_time))
snp_finder.olfh.write('</body>\n</html>\n')
snp_finder.olfh.close()
