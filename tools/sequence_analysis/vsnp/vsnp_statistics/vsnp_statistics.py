#!/usr/bin/env python

import argparse
import gzip
import multiprocessing
import numpy
import os
import pandas
import queue

INPUT_IDXSTATS_DIR = 'input_idxstats'
INPUT_METRICS_DIR = 'input_metrics'
INPUT_READS_DIR = 'input_reads'
OUTPUT_DIR = 'output'
QUALITYKEY = {'!':'0', '"':'1', '#':'2', '$':'3', '%':'4', '&':'5', "'":'6', '(':'7', ')':'8', '*':'9', '+':'10', ',':'11', '-':'12', '.':'13', '/':'14', '0':'15', '1':'16', '2':'17', '3':'18', '4':'19', '5':'20', '6':'21', '7':'22', '8':'23', '9':'24', ':':'25', ';':'26', '<':'27', '=':'28', '>':'29', '?':'30', '@':'31', 'A':'32', 'B':'33', 'C':'34', 'D':'35', 'E':'36', 'F':'37', 'G':'38', 'H':'39', 'I':'40', 'J':'41', 'K':'42', 'L':'43', 'M':'44', 'N':'45', 'O':'46', 'P':'47', 'Q':'48', 'R':'49', 'S':'50', 'T':'51', 'U':'52', 'V':'53', 'W':'54', 'X':'55', 'Y':'56', 'Z':'57', '_':'1', ']':'1', '[':'1', '\\':'1', '\n':'1', '`':'1', 'a':'1', 'b':'1', 'c':'1', 'd':'1', 'e':'1', 'f':'1', 'g':'1', 'h':'1', 'i':'1', 'j':'1', 'k':'1', 'l':'1', 'm':'1', 'n':'1', 'o':'1', 'p':'1', 'q':'1', 'r':'1', 's':'1', 't':'1', 'u':'1', 'v':'1', 'w':'1', 'x':'1', 'y':'1', 'z':'1', ' ':'1'}
READCOLUMNS = ['Sample', 'Reference', 'Fastq File', 'Size', 'Total Reads', 'Mean Read Length', 'Mean Read Quality', 'Reads Passing Q30']
SEP = "\t"


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


def output_read_stats(gzipped, fastq_file, ofh, sampling_number=10000, output_sample=False, dbkey=None, collection=False):
    file_name_base = os.path.basename(fastq_file)
    # Output a 2-column file where column 1 is
    # the labels and column 2 is the values.
    if output_sample:
        # The Sample and Reference columns should be
        # output only once, so this block handles
        # paired reads, where the columns are not
        # output for Read2.
        try:
            # Illumina read file names are something like:
            # 13-1941-6_S4_L001_R1_600000_fastq_gz
            sample = file_name_base.split("_")[0]
        except Exception:
            sample = ""
        # Sample
        ofh.write("Sample%s%s\n" % (SEP, sample))
        ofh.write("Reference%s%s\n" % (SEP, dbkey))
        if collection:
            read = "Read"
        else:
            read = "Read1"
    else:
        read = "Read2"
    # Read
    ofh.write("%s File%s%s\n" % (read, SEP, file_name_base))
    # File Size
    ofh.write("%s File Size%s%s\n" % (read, SEP, nice_size(os.path.getsize(fastq_file))))
    if gzipped.lower() == "true":
        df = pandas.read_csv(gzip.open(fastq_file, "r"), header=None, sep="^")
    else:
        df = pandas.read_csv(open(fastq_file, "r"), header=None, sep="^")
    total_read_count = int(len(df.index) / 4)
    # Readx Total Reads
    ofh.write("%s Total Reads%s%s\n" % (read, SEP, total_read_count))
    # Mean Read Length
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
    ofh.write("%s Mean Read Length%s%s\n" % (read, SEP, "%.1f" % numpy.mean(list_length)))
    # Mean Read Quality
    df_mean = pandas.DataFrame.from_dict(dict_mean, orient='index', columns=['ave'])
    ofh.write("%s Mean Read Quality%s%s\n" % (read, SEP, "%.1f" % df_mean['ave'].mean()))
    # Reads Passing Q30
    reads_gt_q30 = len(df_mean[df_mean['ave'] >= 30])
    reads_passing_q30 = "{:10.2f}".format(reads_gt_q30 / sampling_size)
    ofh.write("%s reads passing Q30%s%s\n" % (read, SEP, reads_passing_q30))
    return total_read_count


def output_statistics(task_queue, read2, collection, gzipped, dbkey, timeout):
    while True:
        try:
            tup = task_queue.get(block=True, timeout=timeout)
        except queue.Empty:
            break
        read_file, idxstats_file, metrics_file, output_file = tup
        total_reads = 0
        with open(output_file, "w") as ofh:
            total_reads += output_read_stats(gzipped, read_file, ofh, output_sample=True, dbkey=dbkey, collection=collection)
            if read2 is not None:
                total_reads += output_read_stats(gzipped, read2, ofh)
            ofh.write("Total Reads%s%d\n" % (SEP, total_reads))
            with open(idxstats_file, "r") as ifh:
                unmapped_reads = 0
                for i, line in enumerate(ifh):
                    items = line.split("\t")
                    if i == 0:
                        # NC_002945.4 4349904 213570 4047
                        ofh.write("All Mapped Reads%s%s\n" % (SEP, items[2]))
                    elif i == 1:
                        # * 0 0 82774
                        unmapped_reads = int(items[3])
                        ofh.write("Unmapped Reads%s%d\n" % (SEP, unmapped_reads))
                percent_str = "Unmapped Reads Percentage of Total"
                if unmapped_reads > 0:
                    unmapped_reads_percentage = "{:10.2f}".format(unmapped_reads / total_reads)
                    ofh.write("%s%s%s\n" % (percent_str, SEP, unmapped_reads_percentage))
                else:
                    ofh.write("%s%s0\n" % (percent_str, SEP))
            with open(metrics_file, "r") as ifh:
                for i, line in enumerate(ifh):
                    if i == 0:
                        # Skip comments.
                        continue
                    items = line.split("\t")
                    if i == 1:
                        # MarkDuplicates 10.338671 98.74%
                        ofh.write("Average Depth of Coverage%s%s\n" % (SEP, items[2]))
                        ofh.write("Reference with Coverage%s%s\n" % (SEP, items[3]))
                    elif i == 2:
                        # VCFfilter 611
                        ofh.write("Good SNP Count%s%s\n" % (SEP, items[1]))
        task_queue.task_done()


def set_num_cpus(num_files, processes):
    num_cpus = int(multiprocessing.cpu_count())
    if num_files < num_cpus and num_files < processes:
        return num_files
    if num_cpus < processes:
        half_cpus = int(num_cpus / 2)
        if num_files < half_cpus:
            return num_files
        return half_cpus
    return processes


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--read1', action='store', dest='read1', required=False, default=None, help='Required: single read')
    parser.add_argument('--read2', action='store', dest='read2', required=False, default=None, help='Optional: paired read')
    parser.add_argument('--dbkey', action='store', dest='dbkey', help='Reference dbkey')
    parser.add_argument('--gzipped', action='store', dest='gzipped', help='Input files are gzipped')
    parser.add_argument('--samtools_idxstats', action='store', dest='samtools_idxstats', required=False, default=None, help='Output of samtools_idxstats')
    parser.add_argument('--output', action='store', dest='output', required=False, default=None, help='Output statisticsfile')
    parser.add_argument('--vsnp_azc', action='store', dest='vsnp_azc', required=False, default=None, help='Output of vsnp_add_zero_coverage')
    parser.add_argument('--processes', action='store', dest='processes', type=int, help='User-selected number of processes to use for job splitting')

    args = parser.parse_args()

    reads_files = []
    idxstats_files = []
    metrics_files = []
    output_files = []
    if args.output is not None:
        # The inputs were not dataset collections, so
        # read1, read2 (possibly) and vsnp_azc will also
        # not be None.
        collection = False
        reads_files.append(args.read1)
        idxstats_files.append(args.samtools_idxstats)
        metrics_files.append(args.vsnp_azc)
        output_files.append(args.output)
    else:
        collection = True
        for file_name in sorted(os.listdir(INPUT_READS_DIR)):
            file_path = os.path.abspath(os.path.join(INPUT_READS_DIR, file_name))
            reads_files.append(file_path)
            base_file_name = get_base_file_name(file_path)
            output_files.append(os.path.abspath(os.path.join(OUTPUT_DIR, "%s.tabular" % base_file_name)))
        for file_name in sorted(os.listdir(INPUT_IDXSTATS_DIR)):
            file_path = os.path.abspath(os.path.join(INPUT_IDXSTATS_DIR, file_name))
            idxstats_files.append(file_path)
        for file_name in sorted(os.listdir(INPUT_METRICS_DIR)):
            file_path = os.path.abspath(os.path.join(INPUT_METRICS_DIR, file_name))
            metrics_files.append(file_path)

    multiprocessing.set_start_method('spawn')
    queue1 = multiprocessing.JoinableQueue()
    num_files = len(output_files)
    cpus = set_num_cpus(num_files, args.processes)
    # Set a timeout for get()s in the queue.
    timeout = 0.05

    for i, output_file in enumerate(output_files):
        read_file = reads_files[i]
        idxstats_file = idxstats_files[i]
        metrics_file = metrics_files[i]
        queue1.put((read_file, idxstats_file, metrics_file, output_file))

    # Complete the output_statistics task.
    processes = [multiprocessing.Process(target=output_statistics, args=(queue1, args.read2, collection, args.gzipped, args.dbkey, timeout, )) for _ in range(cpus)]
    for p in processes:
        p.start()
    for p in processes:
        p.join()
    queue1.join()

    if queue1.empty():
        queue1.close()
        queue1.join_thread()
