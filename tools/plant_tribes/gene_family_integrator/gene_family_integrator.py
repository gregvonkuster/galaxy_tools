#!/usr/bin/env python
import argparse
import subprocess

import utils

OUTPUT_DIR = 'integratedGeneFamilies_dir'

parser = argparse.ArgumentParser()
parser.add_argument('--orthogroup_faa', dest='orthogroup_faa', help="Directory of input fasta datasets")
parser.add_argument('--scaffold', dest='scaffold', help='Orthogroups or gene families proteins scaffold')
parser.add_argument('--method', dest='method', help='Protein clustering method')
parser.add_argument('--orthogroup_fna', dest='orthogroup_fna', default=None, help='Use correspong coding sequences')
parser.add_argument('--output', dest='output', help="Output dataset")
parser.add_argument('--output_dir', dest='output_dir', help="Output dataset file_path directory")

args = parser.parse_args()

# Build the command line.
cmd = 'GeneFamilyIntegrator'
cmd += ' --orthogroup_faa %s' % args.orthogroup_faa
cmd += ' --scaffold %s' % args.scaffold
cmd += ' --method %s' % args.method
if args.orthogroup_fna is not None:
    cmd += ' --orthogroup_fna'
# Run the command.
proc = subprocess.Popen(args=cmd, stderr=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)
rc = proc.wait()
utils.check_execution_errors(rc, proc.stderr)
utils.move_directory_files(OUTPUT_DIR, args.output_dir)
utils.write_html_output(args.output, 'Integrated gene family sequences', args.output_dir)
