#!/usr/bin/env python

import argparse
import gzip
import numpy
import os
import pandas
import shutil

INPUT_IDXSTATS_DIR = 'input_idxstats'
INPUT_METRICS_DIR = 'input_metrics'
INPUT_READS_DIR = 'input_reads'
QUALITYKEY = {'!': '0', '"': '1', '#': '2', '$': '3', '%': '4', '&': '5', "'": '6', '(': '7',
              ')': '8', '*': '9', '+': '10', ',': '11', '-': '12', '.': '13', '/': '14', '0': '15',
              '1': '16', '2': '17', '3': '18', '4': '19', '5': '20', '6': '21', '7': '22',
              '8': '23', '9': '24', ':': '25', ';': '26', '<': '27', '=': '28', '>': '29',
              '?': '30', '@': '31', 'A': '32', 'B': '33', 'C': '34', 'D': '35', 'E': '36',
              'F': '37', 'G': '38', 'H': '39', 'I': '40', 'J': '41', 'K': '42', 'L': '43',
              'M': '44', 'N': '45', 'O': '46', 'P': '47', 'Q': '48', 'R': '49', 'S': '50',
              'T': '51', 'U': '52', 'V': '53', 'W': '54', 'X': '55', 'Y': '56', 'Z': '57',
              '_': '1', ']': '1', '[': '1', '\\': '1', '\n': '1', '`': '1', 'a': '1', 'b': '1',
              'c': '1', 'd': '1', 'e': '1', 'f': '1', 'g': '1', 'h': '1', 'i': '1', 'j': '1',
              'k': '1', 'l': '1', 'm': '1', 'n': '1', 'o': '1', 'p': '1', 'q': '1', 'r': '1',
              's': '1', 't': '1', 'u': '1', 'v': '1', 'w': '1', 'x': '1', 'y': '1', 'z': '1',
              ' ': '1'}


def fastq_to_df(fastq_file, gzipped):
    if gzipped.lower() == "true":
        return pandas.read_csv(gzip.open(fastq_file, "r"), header=None, sep="^")
    else:
        return pandas.read_csv(open(fastq_file, "r"), header=None, sep="^")


def get_base_file_name(file_path):
    base_file_name = os.path.basename(file_path)
    if base_file_name.find(".") > 0:
        # Eliminate the extension.
        return os.path.splitext(base_file_name)[0]
    elif base_file_name.find("_") > 0:
        # The dot extension was likely changed to
        # the " character.
        items = base_file_name.split("_")
        return "_".join(items[0:-1])
    else:
        return base_file_name


def nice_size(size):
    # Returns a readably formatted string with the size
    words = ['bytes', 'KB', 'MB', 'GB', 'TB', 'PB', 'EB']
    prefix = ''
    try:
        size = float(size)
        if size < 0:
            size = abs(size)
            prefix = '-'
    except Exception:
        return '??? bytes'
    for ind, word in enumerate(words):
        step = 1024 ** (ind + 1)
        if step > size:
            size = size / float(1024 ** ind)
            if word == 'bytes':  # No decimals for bytes
                return "%s%d bytes" % (prefix, size)
            return "%s%.1f %s" % (prefix, size, word)
    return '??? bytes'


def output_statistics(reads_files, idxstats_files, metrics_files, output_file, gzipped, dbkey):
    # Produce an Excel spreadsheet that
    # contains a row for each sample.
    columns = ['Reference', 'File Size', 'Mean Read Length', 'Mean Read Quality', 'Reads Passing Q30',
               'Total Reads', 'All Mapped Reads', 'Unmapped Reads', 'Unmapped Reads Percentage of Total',
               'Reference with Coverage', 'Average Depth of Coverage', 'Good SNP Count']
    data_frames = []
    for i, fastq_file in enumerate(reads_files):
        idxstats_file = idxstats_files[i]
        metrics_file = metrics_files[i]
        file_name_base = os.path.basename(fastq_file)
        # Read fastq_file into a data frame.
        fastq_df = fastq_to_df(fastq_file, gzipped)
        total_reads = int(len(fastq_df.index) / 4)
        current_sample_df = pandas.DataFrame(index=[file_name_base], columns=columns)
        # Reference
        current_sample_df.at[file_name_base, 'Reference'] = dbkey
        # File Size
        current_sample_df.at[file_name_base, 'File Size'] = nice_size(os.path.getsize(fastq_file))
        # Mean Read Length
        sampling_size = 10000
        if sampling_size > total_reads:
            sampling_size = total_reads
        fastq_df = fastq_df.iloc[3::4].sample(sampling_size)
        dict_mean = {}
        list_length = []
        for index, row in fastq_df.iterrows():
            base_qualities = []
            for base in list(row.array[0]):
                base_qualities.append(int(QUALITYKEY[base]))
            dict_mean[index] = numpy.mean(base_qualities)
            list_length.append(len(row.array[0]))
        current_sample_df.at[file_name_base, 'Mean Read Length'] = "%.1f" % numpy.mean(list_length)
        # Mean Read Quality
        df_mean = pandas.DataFrame.from_dict(dict_mean, orient='index', columns=['ave'])
        current_sample_df.at[file_name_base, 'Mean Read Quality'] = "%.1f" % df_mean['ave'].mean()
        # Reads Passing Q30
        reads_gt_q30 = len(df_mean[df_mean['ave'] >= 30])
        reads_passing_q30 = "{:10.2f}".format(reads_gt_q30 / sampling_size)
        current_sample_df.at[file_name_base, 'Reads Passing Q30'] = reads_passing_q30
        # Total Reads
        current_sample_df.at[file_name_base, 'Total Reads'] = total_reads
        # All Mapped Reads
        all_mapped_reads, unmapped_reads = process_idxstats_file(idxstats_file)
        current_sample_df.at[file_name_base, 'All Mapped Reads'] = all_mapped_reads
        # Unmapped Reads
        current_sample_df.at[file_name_base, 'Unmapped Reads'] = unmapped_reads
        # Unmapped Reads Percentage of Total
        if unmapped_reads > 0:
            unmapped_reads_percentage = "{:10.2f}".format(unmapped_reads / total_reads)
        else:
            unmapped_reads_percentage = 0
        current_sample_df.at[file_name_base, 'Unmapped Reads Percentage of Total'] = unmapped_reads_percentage
        # Reference with Coverage
        ref_with_coverage, avg_depth_of_coverage, good_snp_count = process_metrics_file(metrics_file)
        current_sample_df.at[file_name_base, 'Reference with Coverage'] = ref_with_coverage
        # Average Depth of Coverage
        current_sample_df.at[file_name_base, 'Average Depth of Coverage'] = avg_depth_of_coverage
        # Good SNP Count
        current_sample_df.at[file_name_base, 'Good SNP Count'] = good_snp_count
        data_frames.append(current_sample_df)
    excel_df = pandas.concat(data_frames)
    excel_file_name = "output.xlsx"
    writer = pandas.ExcelWriter(excel_file_name, engine='xlsxwriter')
    excel_df.to_excel(writer, sheet_name='Sheet1')
    writer.save()
    shutil.move(excel_file_name, output_file)


def process_idxstats_file(idxstats_file):
    all_mapped_reads = 0
    unmapped_reads = 0
    with open(idxstats_file, "r") as fh:
        for i, line in enumerate(fh):
            items = line.split("\t")
            if i == 0:
                # NC_002945.4 4349904 213570 4047
                all_mapped_reads = int(items[2])
            elif i == 1:
                # * 0 0 82774
                unmapped_reads = int(items[3])
    return all_mapped_reads, unmapped_reads


def process_metrics_file(metrics_file):
    ref_with_coverage = '0%'
    avg_depth_of_coverage = 0
    good_snp_count = 0
    with open(metrics_file, "r") as ifh:
        for i, line in enumerate(ifh):
            if i == 0:
                # Skip comments.
                continue
            items = line.split("\t")
            if i == 1:
                # MarkDuplicates 10.338671 98.74%
                ref_with_coverage = items[3]
                avg_depth_of_coverage = items[2]
            elif i == 2:
                # VCFfilter 611
                good_snp_count = items[1]
    return ref_with_coverage, avg_depth_of_coverage, good_snp_count


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--read1', action='store', dest='read1', required=False, default=None, help='Required: single read')
    parser.add_argument('--read2', action='store', dest='read2', required=False, default=None, help='Optional: paired read')
    parser.add_argument('--dbkey', action='store', dest='dbkey', help='Reference dbkey')
    parser.add_argument('--gzipped', action='store', dest='gzipped', help='Input files are gzipped')
    parser.add_argument('--samtools_idxstats', action='store', dest='samtools_idxstats', required=False, default=None, help='Output of samtools_idxstats')
    parser.add_argument('--output', action='store', dest='output', help='Output Excel statistics file')
    parser.add_argument('--vsnp_azc', action='store', dest='vsnp_azc', required=False, default=None, help='Output of vsnp_add_zero_coverage')

    args = parser.parse_args()
    print("args:\n%s\n" % str(args))

    reads_files = []
    idxstats_files = []
    metrics_files = []
    # Accumulate inputs.
    if args.read1 is not None:
        # The inputs are not dataset collections, so
        # read1, read2 (possibly) and vsnp_azc will also
        # not be None.
        reads_files.append(args.read1)
        idxstats_files.append(args.samtools_idxstats)
        metrics_files.append(args.vsnp_azc)
        if args.read2 is not None:
            reads_files.append(args.read2)
            idxstats_files.append(args.samtools_idxstats)
            metrics_files.append(args.vsnp_azc)
    else:
        for file_name in sorted(os.listdir(INPUT_READS_DIR)):
            file_path = os.path.abspath(os.path.join(INPUT_READS_DIR, file_name))
            reads_files.append(file_path)
            base_file_name = get_base_file_name(file_path)
        for file_name in sorted(os.listdir(INPUT_IDXSTATS_DIR)):
            file_path = os.path.abspath(os.path.join(INPUT_IDXSTATS_DIR, file_name))
            idxstats_files.append(file_path)
        for file_name in sorted(os.listdir(INPUT_METRICS_DIR)):
            file_path = os.path.abspath(os.path.join(INPUT_METRICS_DIR, file_name))
            metrics_files.append(file_path)
    output_statistics(reads_files, idxstats_files, metrics_files, args.output, args.gzipped, args.dbkey)
