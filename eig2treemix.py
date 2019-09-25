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
Convert eigensoft *.{ind,geno} files to treemix format.

e.g.
$ ./eig2treemix.py \
        --ind 9908_Bovid749.ind \
        --geno 9908_Bovid749.geno \
        --groups Ovis_aries,Bubalus_bubalis,Bos_grunniens,Bos_Gaurus \
        | gzip -c > 4Bovids.treemix.txt.gz
"""

from __future__ import print_function
import sys
import collections

def parse_geno(filename):
    with open(filename) as f:
        for line in f:
            yield line

def parse_ind(filename):
    lineno = 0
    groups = collections.defaultdict(list)
    with open(filename) as f:
        for line in f:
            fields = line.split()

            #ind = fields[0]
            #sex = fields[1]
            group = fields[2]

            groups[group].append(lineno)
            lineno += 1

    return groups

def parse_args():
    import argparse

    parser = argparse.ArgumentParser(description="Convert eigensoft *.{ind,geno} files to treemix format.")
    parser.add_argument("--ind", required=True, help="eigensoft .ind file")
    parser.add_argument("--geno", required=True, help="eigensoft .geno file")
    parser.add_argument("--groups", help="Comma separated list of groups to include")

    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()

    groups = parse_ind(args.ind)

    if args.groups:
        grouplist = args.groups.split(",")
        for g in grouplist:
            if g not in groups:
                print("group {} not in {}".format(g, args.ind), file=sys.stderr)
    else:
        grouplist = groups.keys()


    # treemix header: space delimited population group names
    print(" ".join(grouplist))

    for line in parse_geno(args.geno):
        out_line = []
        for group in grouplist:
            idx_list =  groups[group]
            a0 = a1 = 0
            for i in idx_list:
                if line[i] == "0":
                    a1 += 2
                elif line[i] == "1":
                    a0 += 1
                    a1 += 1
                elif line[i] == "2":
                    a0 += 2
            out_line.append("{},{}".format(a0, a1))
        print(" ".join(out_line))
