#!/usr/bin/env python

import argparse
import gzip
import numpy
import os
import pandas


QUALITYKEY = {'!':'0', '"':'1', '#':'2', '$':'3', '%':'4', '&':'5', "'":'6', '(':'7', ')':'8', '*':'9', '+':'10', ',':'11', '-':'12', '.':'13', '/':'14', '0':'15', '1':'16', '2':'17', '3':'18', '4':'19', '5':'20', '6':'21', '7':'22', '8':'23', '9':'24', ':':'25', ';':'26', '<':'27', '=':'28', '>':'29', '?':'30', '@':'31', 'A':'32', 'B':'33', 'C':'34', 'D':'35', 'E':'36', 'F':'37', 'G':'38', 'H':'39', 'I':'40', 'J':'41', 'K':'42', 'L':'43', 'M':'44', 'N':'45', 'O':'46', 'P':'47', 'Q':'48', 'R':'49', 'S':'50', 'T':'51', 'U':'52', 'V':'53', 'W':'54', 'X':'55', 'Y':'56', 'Z':'57', '_':'1', ']':'1', '[':'1', '\\':'1', '\n':'1', '`':'1', 'a':'1', 'b':'1', 'c':'1', 'd':'1', 'e':'1', 'f':'1', 'g':'1', 'h':'1', 'i':'1', 'j':'1', 'k':'1', 'l':'1', 'm':'1', 'n':'1', 'o':'1', 'p':'1', 'q':'1', 'r':'1', 's':'1', 't':'1', 'u':'1', 'v':'1', 'w':'1', 'x':'1', 'y':'1', 'z':'1', ' ':'1'}
READCOLUMNS = ['Sample', 'Reference', 'Fastq File', 'Size', 'Total Reads', 'Mean Read Length', 'Mean Read Quality', 'Reads Passing Q30']


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


def output_read_stats(gzipped, fastq_file, fh, sampling_number=10000, output_header=False):
    row_items = []
    try:
        # Illumina read file names are something like:
        # 13-1941-6_S4_L001_R1_600000_fastq_gz
        sample = fastq_file.split("_")[0]
    except Exception:
        sample = ""
    # Sample
    row_items.append(sample)
    # TODO: handle reference instead of setting to n/a.
    row_items.append("n/a")
    # Read
    row_items.append(os.path.basename(fastq_file))
    # File Size
    row_items.append(nice_size(os.path.getsize(fastq_file)))
    header = None
    sep = '^'
    if gzipped:
        df = pandas.read_csv(gzip.open(fastq_file, "r"), header=header, sep=sep)
    else:
        df = pandas.read_csv(open(fastq_file, "r"), header=header, sep=sep)
    total_read_count = int(len(df.index)/4)
    # Fastq Total Reads
    row_items.append(str(total_read_count))
    sampling_size = int(sampling_number)
    if sampling_size > total_read_count:
        sampling_size = total_read_count
    df = df.iloc[3::4].sample(sampling_size)
    dict_mean = {}
    list_length = []
    for index, row in df.iterrows():
        base_qualities = []
        for base in list(row.array[0]):
            base_qualities.append(int(QUALITYKEY[base]))
        dict_mean[index] = numpy.mean(base_qualities)
        list_length.append(len(row.array[0]))
    # Mean Read Length
    row_items.append("%.1f" % numpy.mean(list_length))
    df_mean = pandas.DataFrame.from_dict(dict_mean, orient='index', columns=['ave'])
    # Mean Read Quality
    row_items.append("%.1f" % df_mean['ave'].mean())
    reads_gt_q30 = len(df_mean[df_mean['ave'] >= 30])
    # Reads Passing Q30
    row_items.append(f'{reads_gt_q30/sampling_size:0.1%}')
    # Output read statistics.
    if output_header:
        fh.write("# %s\n" % "\t".join(READCOLUMNS))
    fh.write("%s\n" % "\t".join(row_items))
    return total_read_count


parser = argparse.ArgumentParser()

parser.add_argument('-r1', '--read1', action='store', dest='read1', help='Required: single read')
parser.add_argument('-r2', '--read2', action='store', dest='read2', required=False, default=None, help='Optional: paired read')
parser.add_argument('-gz', '--gzipped', action='store', dest='gzipped', help='Input files are gzipped')
parser.add_argument('-of', '--output', action='store', dest='output', help='Output statisticsfile')

args = parser.parse_args()

total_reads = 0

with open(args.output, "w") as fh:
    total_reads += output_read_stats(args.gzipped, args.read1, fh, output_header=True)
    if args.read2 is not None:
        total_reads += output_read_stats(args.gzipped, args.read2, fh)
