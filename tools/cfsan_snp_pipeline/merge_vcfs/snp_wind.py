#!/usr/bin/env python

import argparse
import os
from os.path import join as j
from itertools import zip_longest


def setup(base_dir, names=[], fwds=[], revs=[], extension='vcf', pattern="{name}.{orient}.{ext}"):
    if fwds and revs and names and len(fwds) != len(revs) != len(names):
        raise ValueError('number of forward reads must equal number of reverse reads and names')
    elif len(fwds) != len(names) or not fwds or not names:
        raise ValueError('number of forward reads must equal number of names')
    with open(j(base_dir, 'snp-unwind.sh'), 'w') as unwind:
        for i, (name, fwd, rev) in enumerate(zip_longest(names, fwds, revs)):
            dir = j(base_dir, str(i))
            sample_dir = j(dir, name)
            os.makedirs(sample_dir)
            target_f = j(sample_dir, pattern.format(name=name, orient=1, ext=extension))
            if rev:
                target_r = j(sample_dir, pattern.format(name=name, orient=2, ext=extension))
            os.symlink(fwd, target_f)
            if rev:
                os.symlink(rev, target_r)
            print(sample_dir)
            if rev:
                unwind.write('unlink {}\n'.format(target_r))
            unwind.write('unlink {}\n'.format(target_f))
            unwind.write('rmdir {}\n'.format(sample_dir))
        unwind.write('rmdir {}\n'.format(dir))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="set up vcf symlink directories for snp-pipeline")
    parser.add_argument('base_dir')
    parser.add_argument('-n', dest='names', type=str, action='append', default=[])
    parser.add_argument('-f', dest='fwds', type=str, action='append', default=[])
    parser.add_argument('-r', dest='revs', type=str, action='append', default=[])
    parser.add_argument('-e', dest='extension', default='vcf')
    parser.add_argument('-p', dest='pattern', default='{name}.{orient}.{ext}')
    params = parser.parse_args()
    setup(**vars(params))
