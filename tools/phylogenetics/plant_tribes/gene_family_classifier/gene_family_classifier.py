#!/usr/bin/env python
import argparse
import os
import shutil

import utils


BUFF_SIZE = 1048576
OUTPUT_DIR = 'geneFamilyClassification_dir'

parser = argparse.ArgumentParser()
parser.add_argument('--input', dest='input', help='Input dataset')
parser.add_argument('--scaffold', dest='scaffold', help='Orthogroups or gene families proteins scaffold')
parser.add_argument('--method', dest='method', help='Protein clustering method')
parser.add_argument('--classifier', dest='classifier', help='Protein classification method')
parser.add_argument('--config_dir', dest='config_dir', help='Directory containing default configuration files')
parser.add_argument('--num_threads', dest='num_threads', type=int, help='Number of threads to use for execution')
parser.add_argument('--super_orthogroups', dest='super_orthogroups', default=None, help='Super orthogroups clustering specification')
parser.add_argument('--single_copy_custom', dest='single_copy_custom', default=None, help='Custom single copy orthogroup configuration')
parser.add_argument('--single_copy_taxa', dest='single_copy_taxa', type=int, default=0, help='Minimum single copy taxa required in orthogroup')
parser.add_argument('--taxa_present', dest='taxa_present', type=int, default=0, help='Minimum taxa required in single copy orthogroup')
parser.add_argument('--orthogroup_fasta', dest='orthogroup_fasta', default=None, help='Flag to create orthogroup sequences')
parser.add_argument('--coding_sequences', dest='coding_sequences', default=None, help='Flag to create orthogroup coding sequences')
parser.add_argument('--save_hmmscan_log', dest='save_hmmscan_log', default=None, help='Flag to save the hmmscan log')
parser.add_argument('--hmmscan_log', dest='hmmscan_log', default=None, help='hmmscan log file')

args = parser.parse_args()

# Build the command line.
cmd = 'GeneFamilyClassifier'
cmd += ' --proteins %s' % args.input
cmd += ' --scaffold %s' % args.scaffold
cmd += ' --method %s' % args.method
cmd += ' --classifier %s' % args.classifier
cmd += ' --config_dir %s' % args.config_dir
cmd += ' --num_threads %d' % args.num_threads
if args.super_orthogroups is not None:
    cmd += ' --super_orthogroups %s' % args.super_orthogroups
if args.single_copy_custom is not None:
    cmd += ' --single_copy_custom %s' % args.single_copy_custom
if args.single_copy_taxa > 0:
    cmd += ' --single_copy_taxa %d' % args.single_copy_taxa
if args.taxa_present > 0:
    cmd += ' --taxa_present %d' % args.taxa_present
if args.orthogroup_fasta is None:
    create_ortho_sequences = False
else:
    create_ortho_sequences = True
    cmd += ' --orthogroup_fasta'
if args.coding_sequences is None:
    create_corresponding_coding_sequences = False
else:
    create_corresponding_coding_sequences = True
    cmd += ' --coding_sequences %s' % args.coding_sequences

# Run the command.
utils.run_command(cmd)

# Handle hmmscan.log output.
if args.classifier in ['hmmscan', 'both']:
    src_hmmscan_log = os.path.join(OUTPUT_DIR, 'hmmscan.log')
    if os.path.exists(src_hmmscan_log):
        if args.save_hmmscan_log is None:
            os.remove(src_hmmscan_log)
        else:
            shutil.move(src_hmmscan_log, args.hmmscan_log)

# Handle orthogroups outputs.
if create_ortho_sequences:
    orthogroups_fasta_src_dir = os.path.join(OUTPUT_DIR, 'orthogroups_fasta')
    orthogroups_fasta_dest_dir = 'output_orthogroups_fasta_dir'
    if not os.path.isdir(orthogroups_fasta_dest_dir):
        os.makedirs(orthogroups_fasta_dest_dir)
    # Remove source direrctory so it won't break dataset collection handler.
    utils.move_directory_files(orthogroups_fasta_src_dir, orthogroups_fasta_dest_dir, remove_source_dir=True)

# Handle single copy orthogroup outputs.
if args.single_copy_custom is not None or args.single_copy_taxa != 0:
    single_copy_fasta_src_dir = os.path.join(OUTPUT_DIR, 'single_copy_fasta')
    single_copy_fasta_dest_dir = 'output_single_copy_fasta_dir'
    if not os.path.isdir(single_copy_fasta_dest_dir):
        os.makedirs(single_copy_fasta_dest_dir)
    # Remove source direrctory so it won't break dataset collection handler.
    utils.move_directory_files(single_copy_fasta_src_dir, single_copy_fasta_dest_dir, remove_source_dir=True)
