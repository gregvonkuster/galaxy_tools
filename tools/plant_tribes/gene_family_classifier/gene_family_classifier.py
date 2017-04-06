#!/usr/bin/env python
import argparse
import os
import shutil
import subprocess
import sys
import tempfile

BUFF_SIZE = 1048576
OUTPUT_DIR = 'geneFamilyClassification_dir'

parser = argparse.ArgumentParser()
parser.add_argument('--input', dest='input', help="Input dataset")
parser.add_argument('--scaffold', dest='scaffold', default='mode', help='Orthogroups or gene families proteins scaffold')
parser.add_argument('--method', dest='method', help='Protein clustering method')
parser.add_argument('--classifier', dest='classifier', help='Protein classification method')
parser.add_argument('--config_dir', dest='config_dir', help='Directory containing default configuration files')
parser.add_argument('--num_threads', dest='num_threads', type=int, help='Number of threads to use for execution')
parser.add_argument('--super_orthogroups', dest='super_orthogroups', default=None, help='Super orthogroups clustering specification')
parser.add_argument('--single_copy_custom', dest='single_copy_custom', default=None, help='Custom single copy orthogroup configuration')
parser.add_argument('--single_copy_taxa', dest='single_copy_taxa', type=int, default=None, help='"Minimum single copy taxa required in orthogroup')
parser.add_argument('--taxa_present', dest='taxa_present', type=int, default=None, help='Minimum taxa required in single copy orthogroup')
parser.add_argument('--orthogroup_fasta', dest='orthogroup_fasta', default=None, help='Flag to create orthogroup sequences')
parser.add_argument('--coding_sequences', dest='coding_sequences', default=None, help='Flag to create orthogroup coding sequences')
parser.add_argument('--save_hmmscan_log', dest='save_hmmscan_log', default=None, help='Flag to save the hmmscan log')
parser.add_argument('--hmmscan_log', dest='hmmscan_log', default=None, help='hmmscan log file')
parser.add_argument('--output_ptortho', dest='output_ptortho', default=None, help='Output for orthogroups')
parser.add_argument('--output_ptortho_dir', dest='output_ptortho_dir', default=None, help='output_ptortho.files_path')
parser.add_argument('--output_ptorthocs', dest='output_ptorthocs', default=None, help='Output for orthogroups with corresponding coding sequences')
parser.add_argument('--output_ptorthocs_dir', dest='output_ptorthocs_dir', default=None, help='output_ptorthocs.files_path')
parser.add_argument('--output_ptsco', dest='output_ptsco', default=None, help='Output for single copy orthogroups')
parser.add_argument('--output_ptsco_dir', dest='output_ptsco_dir', default=None, help='output_ptsco.files_path')

args = parser.parse_args()


def get_stderr_exception(tmp_err, tmp_stderr, tmp_out, tmp_stdout, include_stdout=False):
    tmp_stderr.close()
    # Get stderr, allowing for case where it's very large.
    tmp_stderr = open(tmp_err, 'rb')
    stderr_str = ''
    buffsize = BUFF_SIZE
    try:
        while True:
            stderr_str += tmp_stderr.read(buffsize)
            if not stderr_str or len(stderr_str) % buffsize != 0:
                break
    except OverflowError:
        pass
    tmp_stderr.close()
    if include_stdout:
        tmp_stdout = open(tmp_out, 'rb')
        stdout_str = ''
        buffsize = BUFF_SIZE
        try:
            while True:
                stdout_str += tmp_stdout.read(buffsize)
                if not stdout_str or len(stdout_str) % buffsize != 0:
                    break
        except OverflowError:
            pass
    tmp_stdout.close()
    if include_stdout:
        return 'STDOUT\n%s\n\nSTDERR\n%s\n' % (stdout_str, stderr_str)
    return stderr_str


def move_directory_files(source_dir, destination_dir):
    source_directory = os.path.abspath(source_dir)
    destination_directory = os.path.abspath(destination_dir)
    if not os.path.isdir(destination_directory):
        os.makedirs(destination_directory)
    for dir_entry in os.listdir(source_directory):
        source_entry = os.path.join(source_directory, dir_entry)
        shutil.move(source_entry, destination_directory)


def stop_err(msg):
    sys.stderr.write(msg)
    sys.exit(1)


def write_html_output(output, title, dir):
    with open(output, 'w') as fh:
        fh.write('<html><head><h3>%s</h3></head>\n' % title)
        fh.write('<body><p/><table cellpadding="2">\n')
        fh.write('<tr><th>Size</th><th>Name</th></tr>\n')
        for index, fname in enumerate(sorted(os.listdir(dir))):
            if index % 2 == 0:
                bgcolor = '#D8D8D8'
            else:
                bgcolor = '#FFFFFF'
            try:
                size = str(os.path.getsize(os.path.join(dir, fname)))
            except:
                size = 'unknown'
            link = '<a href="%s" type="text/plain">%s</a>\n' % (fname, fname)
            fh.write('<tr bgcolor="%s"><td>%s</td><td>%s</td></tr>\n' % (bgcolor, size, link))
        fh.write('</table></body></html>\n')


# Define command response buffers.
tmp_out = tempfile.NamedTemporaryFile().name
tmp_stdout = open(tmp_out, 'wb')
tmp_err = tempfile.NamedTemporaryFile().name
tmp_stderr = open(tmp_err, 'wb')
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
if args.single_copy_taxa is not None:
    cmd += ' --single_copy_taxa %d' % args.single_copy_taxa
if args.taxa_present is not None:
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
proc = subprocess.Popen(args=cmd, stderr=tmp_stderr, stdout=tmp_stdout, shell=True)
rc = proc.wait()
# Handle execution errors.
if rc != 0:
    error_message = get_stderr_exception(tmp_err, tmp_stderr, tmp_out, tmp_stdout)
    stop_err( error_message )
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
    if create_corresponding_coding_sequences:
        out_file = args.output_ptorthocs
        orthogroups_fasta_dest_dir = args.output_ptorthocs_dir
        title = 'Orthogroups and corresponding coding sequences files'
    else:
        out_file = args.output_ptortho
        orthogroups_fasta_dest_dir = args.output_ptortho_dir
        title = 'Orthogroups files'
    orthogroups_fasta_src_dir = os.path.join(OUTPUT_DIR, 'orthogroups_fasta')
    move_directory_files(orthogroups_fasta_src_dir, orthogroups_fasta_dest_dir)
    write_html_output(out_file, title, orthogroups_fasta_dest_dir)
# Handle single copy orthogroup outputs.
if args.output_ptsco is not None:
    single_copy_fasta_src_dir = os.path.join(OUTPUT_DIR, 'single_copy_fasta')
    single_copy_fasta_dest_dir = args.output_ptsco_dir
    title = 'Single copy orthogroups files'
    move_directory_files(single_copy_fasta_src_dir, single_copy_fasta_dest_dir)
    write_html_output(args.output_ptsco, title, single_copy_fasta_dest_dir)