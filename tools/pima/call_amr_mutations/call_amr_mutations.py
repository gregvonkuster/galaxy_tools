#!/usr/bin/env python

# NOTE: This tool provides the functionality of both the PIMA filter_varsacn() function
# here https://github.com/appliedbinf/pima_md/blob/main/pima.py#L3012 and the vcf_varscan()
# function here https://github.com/appliedbinf/pima_md/blob/main/pima.py#L3027

import argparse
import os
import subprocess
import sys
import tempfile


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


def filter_varscan(varscan_raw, output):
    cmd = ' '.join(['cat', varscan_raw,
                    '| awk \'(NR > 1 && $9 == 2 && $5 + $6 >= 15)',
                    '{OFS = "\\t";f = $6 / ($5 + $6); gsub(/.*\\//, "", $4);s = $4;gsub(/[+\\-]/, "", s);$7 = sprintf("%.2f%%", f * 100);'
                    'min = 1 / log(length(s) + 2) / log(10) + 2/10;if(f > min){print}}\'',
                    '1>varscan_tmp'])
    run_command(cmd)
    cmd = ' '.join(['cat varscan_tmp',
                    '| awk \'{OFS = "\\t"; print $1,$2,".",$3,$4,-log($14),"PASS",".","GT","1|1"}\'',
                    '1>varscan_vcf'])
    run_command(cmd)
    cmd = ' '.join(['cat varscan_vcf',
                    '| sort -k 1,1 -k 2n,2n',
                    '| awk \'BEGIN{OFS = "\\t";print "##fileformat=VCFv4.2";',
                    'print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE"}{print}\'',
                    '1>' + output])
    run_command(cmd)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--varscan_raw', action='store', dest='varscan_raw', help='Raw varscan mpileup VCF file')
    parser.add_argument('--output', action='store', dest='output', help='Output filtered VCF file')

    args = parser.parse_args()

    filter_varscan(args.varscan_raw, args.output)
