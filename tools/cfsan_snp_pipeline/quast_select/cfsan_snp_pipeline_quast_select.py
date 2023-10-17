#!/usr/bin/env python

import csv
import sys


def pick(rows, key, reverse=False):
    sorted_rows = sorted(rows, key=lambda r: r[key], reverse=reverse)
    return sorted_rows[0]['Assembly']


def int_or_str(token):
    try:
        return int(token)
    except ValueError:
        return str(token)


if __name__ == '__main__':
    path, criterion = sys.argv[1:]
    # QUAST tables have sample info as columns, so we need to transpose the table.
    rows = list(zip(*csv.reader(open(path, "r"), delimiter='\t', dialect='excel')))
    hed = rows.pop(0)
    dict_rows = [{h: int_or_str(r[i]) for i, h in enumerate(hed)} for r in rows]
    if "fewest" in criterion:
        # If it's a count, we want the fewest.
        reverse = False
    else:
        # Otherwise it's a length and we want the longest.
        reverse = True
    print(pick(dict_rows, criterion, reverse))
