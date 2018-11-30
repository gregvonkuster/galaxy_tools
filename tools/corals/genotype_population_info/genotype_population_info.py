#!/usr/bin/env python
"""
Generate a genotype population information file to be used as input to the multilocus_genotype tool.
"""
import argparse
import os
import subprocess

FHEADER = 'header.txt'
FSAMPLE = 'sample.txt'
FSTDERR = 'stderr.txt'
FSTDOUT = 'stdout.txt'

parser = argparse.ArgumentParser()
parser.add_argument('--input_vcf', dest='input_vcf', help='VCF file containing extracted probes from the upstream merged VCF file')
parser.add_argument('--output', dest='output', help='Output dataset'),
args = parser.parse_args()


def check_execution_errors(rc, fstderr, fstdout):
    if rc != 0:
        fh = open(fstdout, 'rb')
        out_msg = fh.read()
        fh.close()
        fh = open(fstderr, 'rb')
        err_msg = fh.read()
        fh.close()
        msg = '%s\n%s\n' % (str(out_msg), str(err_msg))
        stop_err(msg)


def get_response_buffers():
    fstderr = os.path.join(os.getcwd(), FSTDERR)
    fherr = open(fstderr, 'wb')
    fstdout = os.path.join(os.getcwd(), FSTDOUT)
    fhout = open(fstdout, 'wb')
    return fstderr, fherr, fstdout, fhout


def run_command(cmd):
    fstderr, fherr, fstdout, fhout = get_response_buffers()
    proc = subprocess.Popen(args=cmd, stderr=fherr, stdout=fhout, shell=True)
    rc = proc.wait()
    # Check results.
    fherr.close()
    fhout.close()
    check_execution_errors(rc, fstderr, fstdout)


def stop_err(msg):
    sys.exit(msg)


# Extract the header from input_vcf.
cmd = 'grep "#CHROM" %s > %s' % (args.input_vcf, FHEADER)
run_command(cmd)
# Extract the samples based on the header.
cmd = "tr '\t' '\n' < %s > %s" % (FHEADER, FSAMPLE)
run_command(cmd)
# Munge the samples inline.
cmd = "sed -i 1,9d %s" % FSAMPLE
run_command(cmd)
# Generate the output.
cmd = "awk -F'\t' -v OFS='\t' 'NR==0 {print ; next}{print (NR),$0}' %s > %s" % (FSAMPLE, args.output)
run_command(cmd)

