#!/usr/bin/env python

# NOTE: This tool provides the functionality of the PIMA filter_varsacn()
# function here https://github.com/appliedbinf/pima_md/blob/main/pima.py#L3012

import argparse
import re
import subprocess
import sys


def run_command(self, command):
    try:
        return re.split('\\n', subprocess.check_output(command, shell=True).decode('utf-8'))
    except Exception:
        message = 'Command %s failed: exiting...' % command
        sys.exit(message)


def filter_varscan(varscan_raw, output):
    cmd = ' '.join(['cat', varscan_raw,
                    '| awk \'(NR > 1 && $9 == 2 && $5 + $6 >= 15)',
                    '{OFS = "\\t";f = $6 / ($5 + $6); gsub(/.*\\//, "", $4);s = $4;gsub(/[+\\-]/, "", s);$7 = sprintf("%.2f%%", f * 100);'
                    'min = 1 / log(length(s) + 2) / log(10) + 2/10;if(f > min){print}}\'',
                    '1>' + output])
    output = run_command(cmd)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--varscan_raw', action='store', dest='varscan_raw', help='Raw varscan mpileup VCF file')
    parser.add_argument('--output', action='store', dest='output', help='Output filtered VCF file')

    args = parser.parse_args()

    filter_varscan(args.varscan_raw, args.output)
