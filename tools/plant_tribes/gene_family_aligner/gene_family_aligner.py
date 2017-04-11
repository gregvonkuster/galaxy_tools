#!/usr/bin/env python
import argparse
import os
import subprocess

import utils

OUTPUT_DIR = 'geneFamilyAlignments_dir'

parser = argparse.ArgumentParser()
parser.add_argument('--alignment_method', dest='alignment_method', help='Multiple sequence alignments method')
parser.add_argument('--automated_trimming', dest='automated_trimming', default=None, help="Trims alignments using trimAl's ML heuristic trimming approach")
parser.add_argument('--codon_alignments', dest='codon_alignments', default=None, help="Flag for constructing orthogroup multiple codon alignments")
parser.add_argument('--gap_trimming', dest='gap_trimming', default=0, type=float, help='Remove sites in alignments with gaps of')
parser.add_argument('--iterative_realignment', dest='iterative_realignment', type=int, default=0, help='"Maximum number of iterations')
parser.add_argument('--num_threads', dest='num_threads', type=int, help='Number of threads to use for execution')
parser.add_argument('--orthogroup_faa', dest='orthogroup_faa', help="Directory of input fasta datasets")
parser.add_argument('--output', dest='output', help="Output dataset")
parser.add_argument('--output_dir', dest='output_dir', help="Output dataset files_path directory")
parser.add_argument('--pasta_iter_limit', dest='pasta_iter_limit', type=int, default=None, help='"Maximum number of iteration that the PASTA algorithm will execute')
parser.add_argument('--pasta_script_path', dest='pasta_script_path', default=None, help='Path to script for executing pasta')
parser.add_argument('--remove_sequences', dest='remove_sequences', default=0, type=float, help='Remove sequences with gaps of')

args = parser.parse_args()

# Build the command line.
cmd = 'GeneFamilyAligner'
cmd += ' --orthogroup_faa %s' % args.orthogroup_faa
cmd += ' --alignment_method %s' % args.alignment_method
if args.alignment_method == 'pasta':
    if args.pasta_script_path is not None:
        cmd += ' --pasta_script_path %s' % args.pasta_script_path
    if args.pasta_iter_limit is not None:
        cmd += ' --pasta_iter_limit %d' % args.pasta_iter_limit
cmd += ' --num_threads %d' % args.num_threads
if args.codon_alignments is not None:
    cmd += ' --codon_alignments'
if args.automated_trimming is not None:
    cmd += ' --automated_trimming'
if args.gap_trimming > 0:
    cmd += ' --gap_trimming %4f' % args.gap_trimming
if args.remove_sequences > 0:
    cmd += ' --remove_sequences %4f' % args.remove_sequences
if args.iterative_realignment > 0:
    cmd += ' --iterative_realignment %d' % args.iterative_realignment
# Run the command.
proc = subprocess.Popen(args=cmd, stderr=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)
rc = proc.wait()
utils.check_execution_errors(rc, proc.stderr)
if args.codon_alignments is None:
    src_output_dir = OUTPUT_DIR
else:
    src_output_dir = os.path.join(OUTPUT_DIR, 'orthogroups_aln')
utils.move_directory_files(src_output_dir, args.output_dir)
utils.write_html_output(args.output, 'Aligned gene family sequences', args.output_dir)
