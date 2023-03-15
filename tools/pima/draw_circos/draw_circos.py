#!/usr/bin/env python

import argparse
import os
import subprocess
import sys
import tempfile

import Bio.SeqIO
import numpy
import pandas


def format_kmg(number, decimals=0):
    if number == 0:
        return '0'
    magnitude_powers = [10**9, 10**6, 10**3, 1]
    magnitude_units = ['G', 'M', 'K', '']
    for i in range(len(magnitude_units)):
        if number >= magnitude_powers[i]:
            magnitude_power = magnitude_powers[i]
            magnitude_unit = magnitude_units[i]
            return ('{:0.' + str(decimals) + 'f}').format(number / magnitude_power) + magnitude_unit


def load_fasta(fasta_file):
    sequence = pandas.Series(dtype=object)
    for contig in Bio.SeqIO.parse(fasta_file, 'fasta'):
        sequence[contig.id] = contig
    return sequence


def nicenumber(x, round):
    exp = numpy.floor(numpy.log10(x))
    f = x / 10**exp
    if round:
        if f < 1.5:
            nf = 1.0
        elif f < 3.0:
            nf = 2.0
        elif f < 7.0:
            nf = 5.0
        else:
            nf = 10.0
    else:
        if f <= 1.0:
            nf = 1.0
        elif f <= 2.0:
            nf = 2.0
        elif f <= 5.0:
            nf = 5.0
        else:
            nf = 10.0
    return nf * 10.0**exp


def pretty(low, high, n):
    rnge = nicenumber(high - low, False)
    d = nicenumber(rnge / (n - 1), True)
    miny = numpy.floor(low / d) * d
    maxy = numpy.ceil(high / d) * d
    return numpy.arange(miny, maxy + 0.5 * d, d)


def run_command(cmd):
    try:
        tmp_name = tempfile.NamedTemporaryFile(dir=".").name
        tmp_stderr = open(tmp_name, 'wb')
        proc = subprocess.Popen(args=cmd, shell=True, stderr=tmp_stderr.fileno())
        returncode = proc.wait()
        tmp_stderr.close()
        if returncode != 0:
            # Get stderr, allowing for case where it's very large.
            tmp_stderr = open(tmp_name, 'rb')
            stderr = ''
            buffsize = 1048576
            try:
                while True:
                    stderr += tmp_stderr.read(buffsize)
                    if not stderr or len(stderr) % buffsize != 0:
                        break
            except OverflowError:
                pass
            tmp_stderr.close()
            os.remove(tmp_name)
            stop_err(stderr)
    except Exception as e:
        stop_err('Command:\n%s\n\nended with error:\n%s\n\n' % (cmd, str(e)))


def stop_err(msg):
    sys.stderr.write(msg)
    sys.exit(1)


def draw_circos(circos_conf, dnadiff_1coords_file, output_png_dir, reference, reference_sequence_lengths_file, tick_base_conf):
    ofh = open('process_log', 'w')
    ofh.write("circos_conf: %s\n" % str(circos_conf))
    ofh.write("reference: %s\n" % str(reference))
    reference_contigs = reference.index.tolist()
    ofh.write("reference_contigs: %s\n" % str(reference_contigs))
    # Draw one circos plot for each of the contigs in the reference sequence.
    for contig in reference_contigs:
        ofh.write("contig: %s\n" % str(contig))
        contig_dir = os.path.join('circos_dir', contig)
        os.makedirs(contig_dir)
        # Pull the aligned regions out of the dnadiff 1coords file.
        cmd = ' '.join(['cat', dnadiff_1coords_file,
                        '| awk \'$12 == "' + contig + '"\'',
                        '| awk \'{OFS = "\t";print $(NF - 1),$1,$2}\'',
                        '| bedtools complement -g', reference_sequence_lengths_file, '-i -',
                        '| awk \'$3 - $2 >= 25\'',
                        '| bedtools complement -g', reference_sequence_lengths_file, '-i -',
                        '| awk \'{OFS = "\t";print $1,$2,$3}\'',
                        '| awk \'$1 == "' + contig + '"\'',
                        '>', os.path.join(contig_dir, 'alignment.txt')])
        ofh.write("cmd: %s\n" % str(cmd))
        run_command(cmd)
        # Pull the gap regions out of the dnadiff 1coords file.
        cmd = ' '.join(['cat', dnadiff_1coords_file,
                        '| awk \'$12 == "' + contig + '"\'',
                        '| awk \'{OFS = "\t";print $(NF - 1),$1,$2}\'',
                        '| bedtools complement -g', reference_sequence_lengths_file, '-i -',
                        '| awk \'$3 - $2 >= 25\'',
                        '| awk \'{OFS = "\t";print $1,$2,$3}\'',
                        '| awk \'$1 == "' + contig + '"\'',
                        '>', os.path.join(contig_dir, 'gap.txt')])
        ofh.write("cmd: %s\n" % str(cmd))
        run_command(cmd)
        cmd = ' '.join(['cat', reference_sequence_lengths_file,
                        '| awk \'$1 == "' + contig + '"\'',
                        '| awk \'{OFS = "\t";print "chr\t-",$1,$1,0,$2,"plasmid_grey"}\''
                        '>', os.path.join(contig_dir, 'karyotype.txt')])
        ofh.write("cmd: %s\n" % str(cmd))
        run_command(cmd)
        # Figure out the tick labels to use and where to place them.
        # We don't want the last tick since this thing is circular.
        tick_at = pretty(1, len(reference[contig].seq), 12).astype(int)[:-1]
        tick_major = tick_at[1] - tick_at[0]
        tick_minor = tick_major / 5
        cmd = ' '.join(['cat', tick_base_conf,
                        '| awk \'{sub("TICK_MAJOR", "' + str(tick_major) + '", $0);',
                        'sub("TICK_MINOR", "' + str(tick_minor) + '", $0);print}\'',
                        '>', os.path.join(contig_dir, 'tick.conf')])
        ofh.write("cmd: %s\n" % str(cmd))
        run_command(cmd)
        tick_labels = [format_kmg(i) for i in tick_at]
        tick_data = pandas.DataFrame()
        for i in range(len(tick_labels)) :
            tick_data = pandas.concat([tick_data, pandas.Series([contig, tick_at[i], tick_at[i], tick_labels[i]])], axis=1)
        tick_data = tick_data.transpose()
        tick_data.to_csv(path_or_buf=os.path.join(contig_dir, 'tick.txt'), sep='\t', header=False, index=False)
        cmd = ' '.join(['cd', contig_dir,
                        '&& circos --conf', os.path.abspath(circos_conf)])
        ofh.write("cmd: %s\n" % str(cmd))
        run_command(cmd)
        # Move the circos png file for the current
        # contig to the collection output directory.
        ofh.write("\nrenaming: %s to be named %s\n" % (str(os.path.join(contig_dir, 'circos.png')), str(os.path.join(output_png_dir, '%s.png' % contig))))
        os.rename(os.path.join(contig_dir, 'circos.png'), os.path.join(output_png_dir, '%s.png' % contig))
    ofh.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--circos_conf', action='store', dest='circos_conf', help='Circos configuration file')
    parser.add_argument('--dnadiff_1coords_file', action='store', dest='dnadiff_1coords_file', help='Dnadiff 1coords tabular file')
    parser.add_argument('--output_png_dir', action='store', dest='output_png_dir', help='Directory for all circos png outputs')
    parser.add_argument('--reference_file', action='store', dest='reference_file', help='Reference genome fasta file')
    parser.add_argument('--reference_sequence_lengths_file', action='store', dest='reference_sequence_lengths_file', help='Reference sequence lengths tabular file')
    parser.add_argument('--tick_base_conf', action='store', dest='tick_base_conf', help='Tick base configuration file')

    args = parser.parse_args()

    # Load the reference genome into memory.
    reference = load_fasta(args.reference_file)

    draw_circos(args.circos_conf, args.dnadiff_1coords_file, args.output_png_dir, reference, args.reference_sequence_lengths_file, args.tick_base_conf)
