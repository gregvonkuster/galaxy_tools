#!/usr/bin/env python
import argparse
import os

import utils

OUTPUT_DIR = 'geneFamilyPhylogenies_dir'

parser = argparse.ArgumentParser()

parser.add_argument('--alignment_type', dest='alignment_type', help='Input alignments type produced by the GeneFamilyAligner')
parser.add_argument('--bootstrap_replicates', dest='bootstrap_replicates', type=int, default=None, help='Number of replicates for rapid bootstrap analysis')
parser.add_argument('--config_dir', dest='config_dir', help='Directory containing default configuration files')
parser.add_argument('--max_orthogroup_size', dest='max_orthogroup_size', type=int, help='Maximum number of sequences in orthogroup alignments')
parser.add_argument('--method', dest='method', help='Protein clustering method')
parser.add_argument('--min_orthogroup_size', dest='min_orthogroup_size', type=int, help='Minimum number of sequences in orthogroup alignments')
parser.add_argument('--num_threads', dest='num_threads', type=int, help='Number of threads to use for execution')
parser.add_argument('--orthogroup_aln', dest='orthogroup_aln', help='Input dataset files_path')
parser.add_argument('--output_phylip', dest='output_phylip', default=None, help='Output for orthogroup phylip multiple sequence alignments')
parser.add_argument('--output_phylip_dir', dest='output_phylip_dir', default=None, help='output_phylip.files_path')
#parser.add_argument('--output_tree', dest='output_tree', help='Output for phylogenetic trees')
#parser.add_argument('--output_tree_dir', dest='output_tree_dir', help='output_tree.files_path')
parser.add_argument('--rooting_order', dest='rooting_order', default=None, help='Rooting order configuration for rooting trees')
parser.add_argument('--scaffold', dest='scaffold', help='Orthogroups or gene families proteins scaffold')
parser.add_argument('--sequence_type', dest='sequence_type', help='Sequence type used in the phylogenetic inference')
parser.add_argument('--tree_inference', dest='tree_inference', help='Phylogenetic trees inference method')

args = parser.parse_args()

# Build the command line.
cmd = 'GeneFamilyPhylogenyBuilder'
cmd += ' --alignment_type %s' % args.alignment_type
if args.bootstrap_replicates is not None:
    cmd += ' --bootstrap_replicates %d' % args.bootstrap_replicates
cmd += ' --config_dir %s' % args.config_dir
cmd += ' --max_orthogroup_size %d' % args.max_orthogroup_size
cmd += ' --method %s' % args.method
cmd += ' --min_orthogroup_size %d' % args.min_orthogroup_size
cmd += ' --num_threads %d' % args.num_threads
cmd += ' --orthogroup_aln %s' % args.orthogroup_aln
if args.rooting_order is not None:
    cmd += ' --rooting_order %s' % args.rooting_order
cmd += ' --scaffold %s' % args.scaffold
cmd += ' --sequence_type %s' % args.sequence_type
cmd += ' --tree_inference %s' % args.tree_inference

# Run the command.
utils.run_command(cmd)

# Handle outputs.
if args.output_phylip is not None and args.output_phylip_dir is not None:
    src_output_dir = os.path.join(OUTPUT_DIR, 'phylip_aln')
    utils.move_directory_files(src_output_dir, args.output_phylip_dir)
    utils.write_html_output(args.output_phylip, 'Orthogroup phylip multiple sequence alignments', args.output_phylip_dir)
#src_output_dir = os.path.join(OUTPUT_DIR, 'orthogroups_tree')
#utils.move_directory_files(src_output_dir, args.output_tree_dir)
#utils.write_html_output(args.output_tree, 'Phylogenetic trees', args.output_tree_dir)
