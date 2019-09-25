#!/usr/bin/env python
# Copyright (c) 2015 Graham Gower <graham.gower@gmail.com>
#
# Permission to use, copy, modify, and distribute this software for any
# purpose with or without fee is hereby granted, provided that the above
# copyright notice and this permission notice appear in all copies.
#
# THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
# WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
# ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
# WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
# ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
# OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.

"""
Convert eigenstrat format *.{ind,geno} files to fasta-like format.
Output has genotype as 0,1,2, indicating call vs. the reference,
instead of the regular A,C,T,G of a normal fasta file.

e.g.
$ ./eig2fa.py \
        9908_Ancient_Bison.ind \
        9908_Ancient_Bison.geno \
        > 9908_Ancient_Bison.fa
"""

from __future__ import print_function
import re
import sys
import collections

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("usage: {} file.ind file.geno".format(sys.argv[0]), file=sys.stderr)
        sys.exit(1)

    re_space = re.compile(r"[ \t]+")

    individuals = []
    with open(sys.argv[1]) as f:
        for line in f:
            line = line.strip()
            fields = re_space.split(line)
            ind = fields[0]
            sex = fields[1]
            group = fields[2]

            individuals.append("{}_{}".format(ind, group))

    column = collections.defaultdict(list)
    with open(sys.argv[2]) as f:
        for line in f:
            line = line.rstrip("\r\n")
            
            for ind, gt in zip(individuals, line):
                column[ind].append(gt)

    for ind in individuals:
        print(">{}".format(ind))
        gt = "".join(column[ind])
        print(gt)
