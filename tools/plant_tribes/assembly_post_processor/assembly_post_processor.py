#!/usr/bin/env python
import argparse

import utils

parser = argparse.ArgumentParser()
parser.add_argument('--dereplicate', dest='dereplicate', default=None, help='Remove duplicate sequences')
parser.add_argument('--gap_trimming', dest='gap_trimming', type=float, default=0, help='Trim alignments')
parser.add_argument('--gene_family_search', dest='gene_family_search', default=None, help='Targeted gene families')
parser.add_argument('--method', dest='method', default=None, help='Protein clustering method')
parser.add_argument('--min_length', dest='min_length', type=int, default=0, help='Minimum sequence length')
parser.add_argument('--num_threads', dest='num_threads', type=int, help='Number of processors')
parser.add_argument('--prediction_method', dest='prediction_method', help='Coding regions prediction method')
parser.add_argument('--scaffold', dest='scaffold', default=None, help='Gene family scaffold')
parser.add_argument('--score_matrices', dest='score_matrices', default=None, help='Scores matrices')
parser.add_argument('--strand_specific', dest='strand_specific', default=None, help='Strand-specific assembly')
parser.add_argument('--transcripts', dest='transcripts', help='Transcriptome assembly fasta file')

args = parser.parse_args()

# Build the command line.
cmd = 'AssemblyPostProcessor'
if args.dereplicate is not None:
    cmd += ' --dereplicate'
if args.gap_trimming > 0:
    cmd += ' --gap_trimming %4f' % args.gap_trimming
if args.gene_family_search is not None:
    cmd += ' --gene_family_search %s' % args.gene_family_search
if args.method is not None:
    cmd += ' --method %s' % args.method
if args.min_length > 0:
    cmd += ' --min_length %d' % args.min_length
cmd += ' --num_threads %d' % args.num_threads
cmd += ' --prediction_method %s' % args.prediction_method
if args.scaffold is not None:
    cmd += ' --scaffold %s' % args.scaffold
if args.score_matrices is not None:
    cmd += ' --score_matrices %s' % args.score_matrices
if args.strand_specific is not None:
    cmd += ' --strand_specific'
cmd += ' --transcripts %s' % args.transcripts
# Run the command.
utils.run_command(cmd)
